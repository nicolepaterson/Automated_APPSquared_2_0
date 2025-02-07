from Bio.PDB import PDBParser, DSSP
import pandas as pd
import numpy as np
import os
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import BiopythonWarning
import argparse

pd.options.mode.chained_assignment = None  
def parse_args():
    parser = argparse.ArgumentParser(description="Process PDB file and extract residue information.")
    parser.add_argument("-i", "--input_pdb_file", required=True, help="Input PDB file")
    parser.add_argument("-v", "--variant_hash", required=True,
                        help="Alphanumeric character STRING linked to a specific amino acid sequence")
    parser.add_argument("-o", "--out_p", required=True,
                        help="path to output directory")
    return parser.parse_args()

def df_create():
    pdb_filename = args.input_pdb_file
    variant_hash_value = args.variant_hash
    p = PDBParser(QUIET=True)  
    structure = p.get_structure('X', pdb_filename)
    model = structure[0]
    dssp = DSSP(model, pdb_filename)
    data_rows = []
    try:
        res_98_ca_atom = model['A'][98]['CA']
        res_190_ca_atom = model['A'][190]['CA']  
        res_194_ca_atom = model['A'][194]['CA']
        res_153_ca_atom = model['A'][153]['CA']
        distance_to_98_exists = True
        distance_to_190_exists = True
        distance_to_194_exists = True
        distance_to_153_exists = True
    except KeyError:
        distance_to_98_exists = False
        distance_to_190_exists = False
        distance_to_194_exists = False
        distance_to_153_exists = False
    for chain in model:
        for res in chain:
            if not res.id[0].isspace():  # Skip hetero/water residues
                continue
            dssp_key= (chain.id, res.id)
            asa_val, rsa_val= dssp[dssp_key][2], dssp[dssp_key][3] if dssp_key in dssp else (np.nan, np.nan)
            for atom in res:
                distance_to_98_val= (atom.get_vector() - res_98_ca_atom.get_vector()).norm() if distance_to_98_exists else np.nan
                distance_to_190_val= (atom.get_vector() - res_190_ca_atom.get_vector()).norm() if distance_to_190_exists else np.nan
                distance_to_194_val= (atom.get_vector() - res_194_ca_atom.get_vector()).norm() if distance_to_194_exists else np.nan
                distance_to_153_val= (atom.get_vector() - res_153_ca_atom.get_vector()).norm() if distance_to_153_exists else np.nan
                data_rows.append({
                    'Site': res.id[1],
                    'Chain_ID': chain.id,
                    'Node_ID': atom.get_name(),
                    'X_Coord': atom.coord[0],
                    'Y_Coord': atom.coord[1],
                    'Z_Coord': atom.coord[2],
                    'Residue_Type': res.resname,
                    'Distance_to_98' : distance_to_98_val,
                    'Distance_to_190' : distance_to_190_val,
                    'Distance_to_194' : distance_to_194_val,
                    'Distance_to_153' : distance_to_153_val,
                    'RSA' : rsa_val,
                    'ASA' : asa_val,
                    'Variant_Hash': args.variant_hash
                })
            df_residues= pd.DataFrame(data_rows)
        return df_residues
args = parse_args()
#parquet_file_path= f"{args.variant_hash}_epitope.parquet"
csv_file_path= f"{args.out_p}/{args.variant_hash}_epitope.tsv"
df_residues= df_create()
#table_from_df_filtered_residues= pa.Table.from_pandas(df_residues)
#pq.write_table(table_from_df_filtered_residues, parquet_file_path, version='1.0')
print(df_residues)
print(f"Data written to {csv_file_path}")
df_residues.to_csv(csv_file_path, sep="\t", header=False, index=False)

