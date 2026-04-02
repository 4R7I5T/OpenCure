# Clone the CROPSR GitHub repository
git clone https://github.com/iyunsu/cropsr.git

# Navigate to the CROPSR directory
cd cropsr

# Install required Python packages
pip install -r requirements.txt

# Place the sequences.fasta, predictions.gff3 in this directory and run the following:

python cropsr.py -f sequences.fasta -g predictions.gff3 -o output.csv --cas9 -p NGG