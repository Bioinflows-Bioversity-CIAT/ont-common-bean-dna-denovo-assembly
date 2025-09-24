import pandas as pd
from yaml import safe_load
import random
import string

# Snakemake configuration
configfile: "config/config.yaml"

with open(config["resources_config"], "r") as f:
    resources = safe_load(f)

config["organism"]["genome_size_mb"] = int(config["organism"]["genome_size"]/1e6)


sample_units = pd.read_table(config["sample_units"], sep="\t")

base_dir = config["base_dir"]

# Input function,
# return the list of sequencing runs to merge 
# into a single file
def get_read_list_per_sample(wildcards):
    read_file = sample_units[sample_units["sample_name"] == wildcards.sample]["read_file_path"].tolist()
    return read_file

# Input function,
# If exist multiple sequencing runs for one sample
# return the merged file path
def get_raw_reads(wildcards):
    reads = get_read_list_per_sample(wildcards)
    if len(reads) > 1:
        return f"{base_dir}/basecalling/{wildcards.sample}/merged.fastq.gz"
    else:
        return reads[0]

wildcard_constraints:
    sample = "|".join(list(sample_units.sample_name.unique())),
    ref_assebly = "|".join(list(config["scaffolding"]["ref_assembly_name"]))

# After user check the coverage histogram and put the
# thresholds this funciton reads it
def get_read_depth_cutoffs(wildcards):
    depths = pd.read_csv(config['read_depth_path'], sep='\t')
    depths.set_index('sample', drop=False, inplace=True)

    data = depths.loc[wildcards.sample]
    cutoffs = {
        'l' :data.l,
        'm': data.m,
        'h': data.h }
    return cutoffs

def get_big_temp(wildcards):
    # Summary: Get a large temporary directory path for jobs requiring more space.
    # Args: 
    #       wildcards (dict): Snakemake wildcards.
    # Returns:
    #       str: The path to a large temporary directory.
    if config['bigtmp']:
        if config['bigtmp'].endswith("/"):
            return config['bigtmp'] + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
        else:
            return config['bigtmp'] + "/" + "".join(random.choices(string.ascii_uppercase, k=12)) + "/"
    else:
        return tempfile.gettempdir()

def get_prot_ref_assembly(wildcards):
    protein_fasta_path =  config["scaffolding"]["sinteny"][wildcards.ref_assembly]['peptide_seqs']
    return(protein_fasta_path)

def get_gene_annot_ref_assembly(wildcards):
    gene_annot_path = config["scaffolding"]["sinteny"][wildcards.ref_assembly]['gene_annot']
    return(gene_annot_path)

def get_protein_suffix(wildcards):
    suffix = config["scaffolding"]["sinteny"][wildcards.ref_assembly]['prot_sufix']
    return suffix
