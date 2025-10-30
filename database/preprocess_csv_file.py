import pandas as pd
import re
import os
import pickle
import argparse

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--database_csv_file', type=str, default='./database/ICTV_2024.csv')
parser.add_argument('--level', type=str, default='Family')
parser.add_argument('--limited_count', type=int, default=8)
args = parser.parse_args()

def check_uid_format(uid):
    return re.match(r'^[\w]+$', uid) is not None

raw_ictv_data = pd.read_csv(args.database_csv_file)
raw_ictv_data['Family'].fillna(raw_ictv_data['Subfamily'], inplace=True)
raw_ictv_data['Genus'].fillna(raw_ictv_data['Subgenus'], inplace=True)

counts0 = raw_ictv_data[args.level].value_counts()
valid0 = counts0[counts0 >= 0].index
filtered_data0 = raw_ictv_data[raw_ictv_data[args.level].isin(valid0)].copy()
filtered_data0['Virus REFSEQ accession'].fillna(filtered_data0['Virus GENBANK accession'], inplace=True)
filtered_data0 = filtered_data0.dropna(subset=['Virus REFSEQ accession'])

dataset =  set()
ictv_single = []

for index, row in filtered_data0.iterrows():
    if pd.notna(row['Virus REFSEQ accession']):
        dataset.clear() 
        parts = row['Virus REFSEQ accession'].replace(';', ' - ').split(' - ')
        if parts:
            first_part = parts[0]
            if ':' in first_part:
                code = first_part.split(':')[1].strip()
                code = code
            else:
                code = first_part.strip()
        code = re.sub(r'\s*\(.*?\)', '', code)
        if check_uid_format(code):
            row['Virus REFSEQ accession'] = code
            ictv_single.append(row)
    else:
        print(f"The 'Virus REFSEQ accession' value is missing in row {index}.")

ictv_single_df = pd.DataFrame(ictv_single)
ictv_single_df.to_csv('./database/ictv_single.csv', index=False)

def df2txt(df, filename):
    unique_refseq_accessions = set(df['Virus REFSEQ accession'].dropna())
    new_unique_refseq_accessions = set()
    for accession in unique_refseq_accessions:
        parts = accession.replace(';', ' - ').split(' - ')
        for part in parts:
            if part != '':
                if ':' in part:
                    code = part.split(':')[1].strip()
                    new_unique_refseq_accessions.add(code)
                else:
                    new_unique_refseq_accessions.add(part.strip())
    with open(filename, 'w') as file:
        for accession in new_unique_refseq_accessions:
            file.write(accession + '\n')

counts1 = ictv_single_df[args.level].value_counts()
valid1 = counts1[counts1 >= args.limited_count].index
class_counter = {}
taxa_set = set()
filtered_data1 = ictv_single_df[ictv_single_df[args.level].isin(valid1)].copy() 
for index, row in filtered_data1.iterrows():
    taxon = row[args.level]
    taxa_set.add(taxon)
    if class_counter:
        class_id = max(class_counter.values())
    else:
        class_id = 0
    if taxon not in class_counter:
        class_counter[taxon] = class_id + 1
pickle.dump(class_counter, open('./database/class_counter.pkl', "wb" ) )
print(f"count:{len(taxa_set)}")

database_set_df = pd.DataFrame()

for taxon, group_data in filtered_data1.groupby(args.level):
    database_set_df = pd.concat([database_set_df, group_data])

_ = df2txt(database_set_df,'./database/accession.txt')