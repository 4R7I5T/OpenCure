import os
import json
import pandas as pd


def show_file_info(path: str):
    print(f"\n=== {path} ===")
    if not os.path.exists(path):
        print("NOT FOUND")
        return
    size = os.path.getsize(path)
    print(f"Size: {size} bytes")


def preview_csv(path: str, n: int = 5):
    show_file_info(path)
    if not os.path.exists(path):
        return
    try:
        df = pd.read_csv(path)
    except Exception as e:
        print(f"Error reading CSV: {e}")
        return
    print("Columns:", list(df.columns))
    print(df.head(n))


def preview_jsonl(path: str, n: int = 5):
    show_file_info(path)
    if not os.path.exists(path):
        return
    try:
        with open(path, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                if i >= n:
                    break
                line = line.strip()
                if not line:
                    continue
                try:
                    obj = json.loads(line)
                except Exception as e:
                    print(f"Line {i} JSON error: {e}")
                    continue
                print(f"Line {i} keys:", list(obj.keys()))
    except Exception as e:
        print(f"Error reading JSONL: {e}")


def main():
    preview_csv("resp000.csv")
    preview_csv("updated_file.csv")
    preview_csv("output.csv")
    preview_csv("output_annotated.csv")
    preview_csv("cas_offinder_results.tsv")
    preview_jsonl("OpenCureV1.2.jsonl")


if __name__ == "__main__":
    main()
