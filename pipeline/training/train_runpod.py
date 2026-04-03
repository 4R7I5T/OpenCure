#!/usr/bin/env python3
"""
Fine-tune an LLM on OpenCure's CRISPR therapy knowledge using QLoRA.

Designed to run on RunPod (A100 80GB recommended) with Unsloth for 2x
faster training and 60% less VRAM.

Setup (on a RunPod pod with PyTorch template)::

    # Install dependencies
    pip install "unsloth[colab-new] @ git+https://github.com/unslothai/unsloth.git"
    pip install --no-deps "trl<0.9.0" peft accelerate bitsandbytes

    # Set HuggingFace cache to persistent storage
    export HF_HOME=/runpod-volume/huggingface

    # Generate training data (from the OpenCure repo)
    python -m pipeline.training.export_training_data -o /runpod-volume/data/opencure_train.jsonl

    # Run training
    python pipeline/training/train_runpod.py \\
        --data /runpod-volume/data/opencure_train.jsonl \\
        --output /runpod-volume/outputs/opencure-bio-llm

    # Or with a different base model:
    python pipeline/training/train_runpod.py \\
        --model Qwen/Qwen2.5-7B-Instruct \\
        --data /runpod-volume/data/opencure_train.jsonl

Estimated cost: $3–10 on RunPod (A100 80GB, ~45–90 min).
"""

import argparse
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Fine-tune an LLM on OpenCure CRISPR therapy data (QLoRA)",
    )
    parser.add_argument(
        "--model", default="meta-llama/Meta-Llama-3.1-8B-Instruct",
        help="HuggingFace model ID (default: LLaMA 3.1 8B Instruct)",
    )
    parser.add_argument(
        "--data", required=True,
        help="Path to training data JSONL (Alpaca format)",
    )
    parser.add_argument(
        "--output", "-o", default="./opencure-bio-llm",
        help="Output directory for adapter weights",
    )
    parser.add_argument(
        "--epochs", type=int, default=3,
        help="Number of training epochs (default: 3)",
    )
    parser.add_argument(
        "--lr", type=float, default=2e-4,
        help="Learning rate (default: 2e-4)",
    )
    parser.add_argument(
        "--batch-size", type=int, default=4,
        help="Per-device batch size (default: 4)",
    )
    parser.add_argument(
        "--grad-accum", type=int, default=4,
        help="Gradient accumulation steps (default: 4, effective batch=16)",
    )
    parser.add_argument(
        "--lora-r", type=int, default=64,
        help="LoRA rank (default: 64)",
    )
    parser.add_argument(
        "--lora-alpha", type=int, default=128,
        help="LoRA alpha (default: 128)",
    )
    parser.add_argument(
        "--max-seq-len", type=int, default=2048,
        help="Maximum sequence length (default: 2048)",
    )
    parser.add_argument(
        "--merge", action="store_true",
        help="Merge LoRA weights into base model after training",
    )
    parser.add_argument(
        "--val-split", type=float, default=0.05,
        help="Fraction of data for validation (default: 0.05)",
    )
    args = parser.parse_args()

    # ---- Imports (heavy, so deferred) ----
    try:
        from unsloth import FastLanguageModel
    except ImportError:
        sys.exit(
            "Unsloth not installed. Run:\n"
            '  pip install "unsloth[colab-new] @ git+https://github.com/unslothai/unsloth.git"\n'
            '  pip install --no-deps "trl<0.9.0" peft accelerate bitsandbytes'
        )

    import torch
    from trl import SFTTrainer
    from transformers import TrainingArguments
    from datasets import load_dataset

    print(f"[OpenCure] Loading base model: {args.model}")
    print(f"[OpenCure] Training data: {args.data}")
    print(f"[OpenCure] Output: {args.output}")
    print(f"[OpenCure] GPU: {torch.cuda.get_device_name(0)}")
    print(f"[OpenCure] VRAM: {torch.cuda.get_device_properties(0).total_mem / 1e9:.1f} GB")
    print()

    # ---- 1. Load model with 4-bit quantization ----
    model, tokenizer = FastLanguageModel.from_pretrained(
        model_name=args.model,
        max_seq_length=args.max_seq_len,
        dtype=None,  # auto-detect (bf16 on Ampere+)
        load_in_4bit=True,
        device_map="auto",
    )

    # ---- 2. Add LoRA adapters ----
    model = FastLanguageModel.get_peft_model(
        model,
        r=args.lora_r,
        target_modules=[
            "q_proj", "k_proj", "v_proj", "o_proj",
            "gate_proj", "up_proj", "down_proj",
        ],
        lora_alpha=args.lora_alpha,
        lora_dropout=0.05,
        bias="none",
        use_gradient_checkpointing="unsloth",
    )

    trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    total = sum(p.numel() for p in model.parameters())
    print(f"[OpenCure] Trainable parameters: {trainable:,} / {total:,} "
          f"({trainable/total*100:.2f}%)")

    # ---- 3. Load and format dataset ----
    dataset = load_dataset("json", data_files=args.data, split="train")
    print(f"[OpenCure] Dataset: {len(dataset)} examples")

    # Alpaca prompt template (matches LLaMA 3 instruction format)
    TEMPLATE = (
        "Below is an instruction that describes a task, paired with further "
        "context. Write a response that appropriately completes the request.\n\n"
        "### Instruction:\n{instruction}\n\n"
        "### Input:\n{input}\n\n"
        "### Response:\n{output}"
    )

    def format_prompts(examples):
        texts = []
        for inst, inp, out in zip(
            examples["instruction"], examples["input"], examples["output"]
        ):
            text = TEMPLATE.format(instruction=inst, input=inp or "", output=out)
            texts.append(text + tokenizer.eos_token)
        return {"text": texts}

    dataset = dataset.map(format_prompts, batched=True)

    # Train/val split
    if args.val_split > 0:
        split = dataset.train_test_split(test_size=args.val_split, seed=42)
        train_dataset = split["train"]
        eval_dataset = split["test"]
        print(f"[OpenCure] Train: {len(train_dataset)}, Val: {len(eval_dataset)}")
    else:
        train_dataset = dataset
        eval_dataset = None

    # ---- 4. Training ----
    os.makedirs(args.output, exist_ok=True)

    training_args = TrainingArguments(
        per_device_train_batch_size=args.batch_size,
        gradient_accumulation_steps=args.grad_accum,
        warmup_ratio=0.03,
        num_train_epochs=args.epochs,
        learning_rate=args.lr,
        fp16=not torch.cuda.is_bf16_supported(),
        bf16=torch.cuda.is_bf16_supported(),
        logging_steps=10,
        optim="adamw_8bit",
        weight_decay=0.01,
        lr_scheduler_type="cosine",
        output_dir=args.output,
        save_strategy="epoch",
        eval_strategy="epoch" if eval_dataset else "no",
        load_best_model_at_end=bool(eval_dataset),
        report_to="none",
        max_grad_norm=1.0,
    )

    trainer = SFTTrainer(
        model=model,
        tokenizer=tokenizer,
        train_dataset=train_dataset,
        eval_dataset=eval_dataset,
        dataset_text_field="text",
        max_seq_length=args.max_seq_len,
        packing=True,
        args=training_args,
    )

    print("[OpenCure] Starting training...")
    result = trainer.train()
    print(f"[OpenCure] Training complete. Loss: {result.training_loss:.4f}")

    # ---- 5. Save ----
    adapter_path = os.path.join(args.output, "adapter")
    model.save_pretrained(adapter_path)
    tokenizer.save_pretrained(adapter_path)
    print(f"[OpenCure] LoRA adapter saved to {adapter_path}")

    if args.merge:
        merged_path = os.path.join(args.output, "merged")
        print(f"[OpenCure] Merging weights into full model at {merged_path}...")
        model.save_pretrained_merged(merged_path, tokenizer, save_method="merged_16bit")
        print(f"[OpenCure] Merged model saved to {merged_path}")

    # ---- 6. Save training metadata ----
    meta = {
        "base_model": args.model,
        "lora_r": args.lora_r,
        "lora_alpha": args.lora_alpha,
        "epochs": args.epochs,
        "learning_rate": args.lr,
        "effective_batch_size": args.batch_size * args.grad_accum,
        "train_examples": len(train_dataset),
        "val_examples": len(eval_dataset) if eval_dataset else 0,
        "trainable_params": trainable,
        "total_params": total,
        "final_loss": result.training_loss,
        "gpu": torch.cuda.get_device_name(0),
    }
    with open(os.path.join(args.output, "training_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print("[OpenCure] Done.")


if __name__ == "__main__":
    main()
