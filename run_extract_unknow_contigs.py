import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'contigs.fasta')
args = parser.parse_args()

pred_df = pd.read_csv("./result/prediction.csv")
existing_contigs = set(pred_df["contig_name"])
output_records = []
for record in SeqIO.parse(args.contigs):
    if record.id not in existing_contigs:
        output_records.append(record)
SeqIO.write(output_records, "result/unknow_contigs.fasta", "fasta")
print(f"Extraction completed! {len(output_records)} sequences have been saved to result/unknown_contigs.fasta.")