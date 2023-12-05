configfile: 'config_appsquared.yaml'

from datetime import date
import os
import re

rule all:
    input: 
        pdb = expand("/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_pdbs/{today}/{variant_hash}.pdb",today=config['today'],variant_hash=config['vhash_list'])
rule simple:
    input: 
        fasta = "/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_fasta/{today}/{variant_hash}.fasta"
    output: 
        pdb = "/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_pdbs/{today}/{variant_hash}.pdb"
    threads:
        16 
    shell: """python3 /scicomp/groups/OID/NCIRD/ID/VSDB/GAT/docker_containers/alphafold/alphafold-2.2.4/docker/run_docker.py \
--fasta_paths={input.fasta} \
--model_preset='multimer' \
--max_template_date='2020-05-14' \
--data_dir='/scicomp/reference/alphafold' \
--output_dir="/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_pdbs/{config[today]}"
""" 
#--output_dir='{/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_pdbs/{config[today]}}'
#""" 
