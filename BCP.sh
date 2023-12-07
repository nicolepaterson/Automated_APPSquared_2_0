#!usr/bin/bash
# @nicolemariepaterson
# HA pipeline

set -e
set -u
set -o pipefail

timestamp=$(date +%Y-%m-%d)
echo $timestamp

#paths not expected to change frequently for input/output
gwa="/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/"
gwa_output="/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_archive/$timestamp"
pdb_dir="/scicomp/groups/OID/NCIRD/ID/VSDB/GAT/cdp_pdbs/2023-12-04"
home="/scicomp/home-pure/qxa4/BCP_V2"

#run prepwizard if not already run on structures, optimizes sidechain conformations and finds energy minima for structure
function run_prepwizard { 
ml schrodinger/2020-2
filename=$(basename $file)
echo $filename
cd $pdb_dir
run utilities/prepwizard -disulfides $filename prep_$filename
cd $home
}

#run Rosetta energy score for ddG, prep file for cdp upload
function run_rosetta { 
source /scicomp/home-pure/qxa4/miniconda3/etc/profile.d/conda.sh
source /scicomp/home-pure/qxa4/miniconda3/bin/activate /scicomp/home-pure/qxa4/.conda/envs/appsquared
filename=$(basename "$file" .pdb)
echo $filename
module load rosetta/3.13-cluster
timestamp=$(date +%Y-%m-%d)
mkdir $gwa_output
score_jd2.default.linuxgccrelease -in:file:s $file -out:file:scorefile $gwa_output/$filename.scorefile.txt
sed -i "s/SCORE: //g" $gwa_output/$filename.scorefile.txt
python dataframe.py -i $gwa_output/$filename.scorefile.txt -v $filename -t $timestamp
sed -i 's/"//g' $gwa_output/$filename.scorefile.csv
sed -n '1p' $gwa_output/$filename.scorefile.csv
}

#computes the glycosylation distances for each site from pdb file
function run_glyc_dist { 
filename=$(basename $file)
variant_hash=$(basename $file .pdb)
source /scicomp/home-pure/qxa4/miniconda3/bin/conda
source /scicomp/home-pure/qxa4/miniconda3/bin/activate --stack /scicomp/home-pure/qxa4/.conda/envs/glyc
bash run_glyc.sh -i $pdb_dir -t $timestamp -o $gwa_output/$filename.ASA_glyc.csv
}

#runs getcontacts and compiles the atomic contact types, formats for cdp upload
#source ~/miniconda3/etc/profile.d/conda.sh
#source /scicomp/home-pure/qxa4/miniconda3/bin/conda
function run_get_contacts {
source /scicomp/home-pure/qxa4/miniconda3/bin/activate --stack /scicomp/home-pure/qxa4/.conda/envs/getcontacts
filename=$(basename $file)
variant_hash=$(basename $file .pdb)
echo $variant_hash
python get_static_contacts.py --structure $pdb_dir/$filename --output $gwa_output/$variant_hash.static_contacts.csv --itypes all
}

function run_get_contacts_clean {
for file in $gwa_output/*static_contacts.csv
do
    #filename=$(basename "$file")
    filename=$(basename "$file")
    variant_hash=$(basename $file .static_contacts.csv)
    python getcontacts_df.py -i $file -v $variant_hash -t $timestamp
#    filename=$(basename "$file")
done
}

#generates a .mae file and optimizes for docking runs
function convert_pdb_mae {
run pdbconvert -ipdb $gwa_output/$filename -omae $gwa_output/$filename.mae
}

function prep_mae {
run /utilities/prepwizard $gwa_output/$filename.mae prep_$filename.mae
}

#Generates a control file for grid generation for Glide ligand dock
function generate_gridfile {
filename=$(basename "$file" .mae);
gridtable="FORCEFIELD   OPLS_2005\nGRID_CENTER_ASL residue.num 153\nGRIDFILE    $gwa_output/"$filename"grid.zip\nINNERBOX   10, 10, 10\nOUTERBOX   30, 30, 30\nRECEP_FILE   "$file""
echo -e $gridtable > $gwa_output/$filename"grid.inp"
}
#Generate a control file for running Glide sialic acid docks from grid file
function run_grid_gen {
filename=$(basename "$file" .mae);
docktable="FORCEFIELD   OPLS_2005\nGRIDFILE   $gwa_output/"$filename"grid.zip\nLIGANDFILE   /scicomp/home-pure/qxa4/BCP_V2/ligprep_6-sialyl-out.maegz\nNREPORT   2\nPOSTDOCK_XP_DELE   0.5\nPRECISION   XP\nWRITE_XP_DESC   True\nWRITE_RES_INTERACTION   True\nWRITE_CSV       True"
echo -e $docktable > $gwa_output/$filename.dock.inp
}

#Runs the grid gen and sialic acid docking from the control file
function run_sialic_acid_docking {
filename=$(basename "$file" .mae);
run glide $gwa_output/$filename -OVERWRITE -HOST localhost -TMPLAUNCHDIR -LOCAL
}

#Generate a control file for antibody docks and run antibody dock
function run_antibody_docking {
run -FROM psp piper.py -jobname $filename_dock -receptor prep_7TZ5.pdb -receptor_chain B,C -poses 10 -rotations 70000 -OMPI 1 -JOBID -antibody -ligand $filename -ligand_chain A -use_nonstandard_residue yes -HOST localhost:1 -TMPLAUNCHDIR
}

#Runs the pipeline in order
for file in $pdb_dir/*.pdb;
do
echo $timestamp "RUNNING PREPWIZARD"
run_prepwizard 
done

for file in $pdb_dir/prep*.pdb
do
echo $timestamp "RUNNING ROSETTA"
run_rosetta 
done

for file in $pdb_dir/prep*.pdb
do
echo $timestamp "RUNNING GLYCOSYLATION DISTANCE CALCULATOR"
run_glyc_dist 
done

for file in $pdb_dir/prep*.pdb
do
echo $timestamp "CALCULATING INTERMOLECULAR CONTACTS"
run_get_contacts | run_get_contacts_clean 
done

for file in $pdb_dir/prep*.pdb
do
echo $timestamp "DOCKING TO SIALIC ACID"
convert_pdb_mae | prep_mae | generate_gridfile | run_grid_gen | run_sialic_acid_docking |
done

#echo $timestamp "SUBMITTING ANTIBODY DOCKING"
#run_antibody_docking
