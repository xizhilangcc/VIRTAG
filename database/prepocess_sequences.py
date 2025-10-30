import pandas as pd
from Bio import SeqIO
import pickle
from shutil import which
import subprocess
import argparse


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--level', type=str, default='Family')
args = parser.parse_args()

folder_path = './database'
data = pd.read_csv(f'./database/ictv_single.csv')
fasta_file = f'./database/database.fasta'
output_records = []
class_counter = pickle.load(open('./database/class_counter.pkl', 'rb'))


for record in SeqIO.parse(fasta_file, 'fasta'):
    seq_name = record.id.rsplit('.', 1)[0]
    matching_row = data[data['Virus REFSEQ accession'] == seq_name]
    if not matching_row.empty:
        family = matching_row['Family'].values[0]
        genus = matching_row['Genus'].values[0]
    if not matching_row.empty:
        if args.level == 'Family':
            family_name = matching_row['Family'].values[0]
            class_number = class_counter.get(family_name, None)
            if class_number is not None:
                new_id = f"{family_name}_family_{class_number}_{seq_name}"
        elif args.level == 'Genus':
            genus_name = matching_row['Genus'].values[0]
            class_number = class_counter.get(genus_name, None)  # 继续按属计数
            if class_number is not None:
                new_id = f"{genus_name}_genus_{class_number}_{seq_name}"
        elif args.level == 'Species':
            species_name = matching_row['Species'].values[0]
            class_number = class_counter.get(species_name, None)
            if class_number is not None:
                new_id = f"{species_name}_species_{class_number}_{seq_name}"
        else:
            continue
        record.id = new_id
        record.description = ""
        output_records.append(record)


with open(f'{folder_path}/database_with_class.fasta', 'w') as out_file:
    SeqIO.write(output_records, out_file, 'fasta')

threads = 8
prodigal = "prodigal"
if which("pprodigal") is not None:
    print("Using parallelized prodigal...")
    prodigal = f'pprodigal -T {threads}'
prodigal_cmd = f'{prodigal} -i {folder_path}/database_with_class.fasta -a {folder_path}/database_with_class_prodigal_protein.fasta -f gff -p meta'
print("Running prodigal...")
_ = subprocess.check_call(prodigal_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
fasta_file = f"{folder_path}/database_with_class_prodigal_protein.fasta"
protein_ids = set()
for record in SeqIO.parse(fasta_file, "fasta"):
    protein_ids.add(record.id)
combined_list = list(protein_ids)
protein_id = sorted(combined_list)
contig_id = [item.rsplit("_", 1)[0] for item in protein_id]
description = [item.replace(".", "").replace("_", "") for item in protein_id]
gene2genome = pd.DataFrame({"protein_id": protein_id, "contig_id": contig_id ,"keywords": description})
gene2genome.to_csv(f"{folder_path}/ALL_gene_to_genomes.csv", index=None)