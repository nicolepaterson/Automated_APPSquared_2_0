# automated APPSquared
This version of the pipeline is intended to be run on a weekly schedule. The pipeline first pulls variant hashes with counts >10 from the global variant hash alignments tables, generates an amino acid fasta file from these variant hashes (if it fails to generate a sequence, aa_seq_fails.txt records this information in the local directory). The variant_hash_tracking file in this repo is used for tracking purposes as the name implies, and this file is pulled, updated, and pushed back to the repo after a run. The protein_modeling.isolate_name table is updated with the isolate name that corresponds to gaa.isolate_ranking=1 and updated in CDP for consistency in naming conventions. The amino acid fastas are run in AlphaFold2 using a snakemake workflow with variable control resource allocation in the scbs-vsdb-01 server.  It then runs the APPSquared pipeline on these structures (https://pypi.org/project/appsquared/) archives pdbs and raw output in the GAT group working area and securely copies formatted data tables for upload to CDP. 

![appsquared_workflow]([appsquared_workflow.png](https://github.com/nicolepaterson/Automated_APPSquared_2_0/blob/main/appsquared_workflow.png))

To install the conda environments:
```
conda env create --name appsquared --file=appsquared.yaml
conda env create --name glyc --file=glyc.yaml
conda env create --name getcontacts --file=getcontacts.yaml
```

To run the pipeline for surveillance purposes:
```
bash run_docker_AF.sh -o </path/to/output/dir> -s <start_date format: %Y-%m-%d>
```

        example: bash run_docker_AF.sh -o /home/nicole -s 2023-10-01


pipeline usage standalone (does not generate AlphaFold2 structures or upload to CDP):
```
bash BCP.sh -d </path/to/pdb/files> -n name_of_output_dir
```
table uploads, run from cdp-client-02:
```
./SIMPLE_auto do
```

## Authors and acknowledgment
Nicholas Kovacs, Brian Mann, Kristine Lacek, Norman Hassell, Matthew Wersebe, Sam Shepard have all made meaningful contributions to this project in the form of contributing code that was either used directly as noted under the shebang in each file or modified for this purpose.

## License
open source

## Project status
Beta
# Automated_APPSquared_2_0
Automated version of APPSquared Pipeline plus updates
