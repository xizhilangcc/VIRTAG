import numpy as np
import pandas as pd
import os
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import subprocess
import argparse
import re
VERSION = "VIRTAG version : 1.0"
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--contigs', type=str, default = 'result/unknow_contigs.fasta')
parser.add_argument('--len', type=int, default=1700)
parser.add_argument('--outpath', type=str, default= "result_unknow")
parser.add_argument('--database_path', type=str, default= "./database")
parser.add_argument('--class_counter_path', type=str, default = "./database/class_counter.pkl")
parser.add_argument('-v', '--version', action='version', version=VERSION)
args = parser.parse_args()

folder_path = args.database_path
if os.path.exists(f"{args.outpath}"):
    print("folder {0} exist... please make sure the result folder is different when you run it multiple times".format(f"{args.outpath}"))

def check_folder(file_name):
    if not os.path.exists(file_name):
        _ = os.makedirs(file_name)
    else:
        print("folder {0} exist... cleaning dictionary".format(file_name))
        if os.listdir(file_name):
            try:
                _ = subprocess.check_call("rm -rf {0}".format(file_name), shell=True)
                _ = os.makedirs(file_name)
                print("Dictionary cleaned")
            except:
                print("Cannot clean your folder... permission denied")
                exit(1)
check_folder(f"{args.outpath}")
check_folder(f"{args.outpath}/input")
check_folder(f"{args.outpath}/pred")
check_folder(f"{args.outpath}/Split_files")
check_folder(f"{args.outpath}/network")
check_folder(f"{args.outpath}/Cyber_data")

if os.path.exists(f"{folder_path}/database.self-diamond.tab.abc"):
    if os.path.getsize(f"{folder_path}/database.self-diamond.tab.abc") == 0:
        try:
            make_diamond_cmd = f'diamond makedb --threads 128 --in {folder_path}/database_with_class_prodigal_protein.fasta -d {folder_path}/database.dmnd'
            print("Creating Diamond database...")
            _ = subprocess.check_call(make_diamond_cmd, shell=True)
            diamond_cmd = 'diamond blastp --threads 128 --sensitive -d {folder_path}/database.dmnd -q {folder_path}/database_with_class_prodigal_protein.fasta -o {folder_path}/database.self-diamond.tab'
            print("Running Diamond...")
            _ = subprocess.check_call(diamond_cmd, shell=True)
            diamond_out_fp = f"{folder_path}/database.self-diamond.tab"
            database_abc_fp = f"{folder_path}/database.self-diamond.tab.abc"
            _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
        except:
            print("create database failed")
            exit(1)
else:
    try:
        make_diamond_cmd = f'diamond makedb --threads 128 --in {folder_path}/database_with_class_prodigal_protein.fasta -d {folder_path}/database.dmnd'
        print("Creating Diamond database...")
        _ = subprocess.check_call(make_diamond_cmd, shell=True)

        diamond_cmd = f'diamond blastp --threads 128 --sensitive -d {folder_path}/database.dmnd -q {folder_path}/database_with_class_prodigal_protein.fasta -o {folder_path}/database.self-diamond.tab'
        print("Running Diamond...")
        _ = subprocess.check_call(diamond_cmd, shell=True)
        diamond_out_fp = f"{folder_path}/database.self-diamond.tab"
        database_abc_fp = f"{folder_path}/database.self-diamond.tab.abc"
        _ = subprocess.check_call("awk '$1!=$2 {{print $1,$2,$11}}' {0} > {1}".format(diamond_out_fp, database_abc_fp), shell=True)
    except:
        print("create database failed")
        exit(1)

#####################################################################
##########################    Start Program  ########################
#####################################################################


def special_match(strg):
    cleaned_str = re.sub(r'[^ATGCN]', '', strg)
    return cleaned_str


cnt = 0
file_id = 0
records = []
for record in SeqIO.parse(args.contigs, 'fasta'):
    if cnt !=0 and cnt%300000 == 0:
        SeqIO.write(records, f"{args.outpath}/Split_files/contig_"+str(file_id)+".fasta","fasta") 
        records = []
        file_id+=1
        cnt = 0
    seq = str(record.seq)
    record.seq = Seq(special_match(seq.upper()))
    if len(record.seq) > args.len:
        records.append(record)
        cnt+=1

SeqIO.write(records, f"{args.outpath}/Split_files/contig_"+str(file_id)+".fasta","fasta")
file_id+=1


for i in range(file_id):
    cmd = f"mv {args.outpath}/Split_files/contig_"+str(i)+f".fasta {args.outpath}/input/"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Moving file Error for file {0}".format("contig_"+str(i)))
        continue

    cmd = f"python 'run_construct_graph_network_unknow.py' --database_path {args.database_path} --outpath {args.outpath} --n " + str(i)
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("Knowledge Graph Error for file {0}".format("contig_"+str(i)))
        cmd = f"rm {args.outpath}/input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = f"python 'run_cluster.py' --database_path {args.database_path} --class_counter_path {args.class_counter_path} --outpath {args.outpath}"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("label propagation Error for file {0}".format("contig_"+str(i)))
        cmd = f"rm {args.outpath}/input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue

    cmd = f"rm {args.outpath}/input/*"
    try:
        out = subprocess.check_call(cmd, shell=True)
    except:
        print("rm Error for file {0}".format("contig_"+str(i)))
        cmd = f"rm {args.outpath}/input/*"
        out = subprocess.check_call(cmd, shell=True)
        continue