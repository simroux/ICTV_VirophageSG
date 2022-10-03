# ICTV_VirophageSG
Data and information related to work by the ICTV Virophage Study Group

# Virophage affiliation scripts
Set of scripts to search for virophage marker genes and assign new sequences to existing virophage taxa

## Dependencies
These scripts rely on prodigal, hmmer, and blast. We recommend installing these dependencies in a conda environment with mamba, as follows:
```
mamba create -n virophage_classification -c bioconda prodigal hmmer blast perl 
```
or
```
mamba env create -f virophage_classification.yaml
```

## Running the script
The script requires two arguments: an input file of the genome (fasta), and a directory for the output files.
You can run it on the test sequences (obtained from Bellas & Sommaruga, 2021, https://doi.org/10.1186/s40168-020-00956-0), as follows:
```
mkdir Test_output_new
./run_virophage_affi.pl -i Test_input/Input_test_Bellas-Sommaruga_GKS.fna -o Test_output_new/
```

This should creates files similar to the output currently provided in "Test_output"

## Output

The main output file is XXXX_final_affiliation.tsv, which assigns each input sequence to a virophage taxon, if possible.

The file XXXX_best_markers-step_1.tsv lists the specific predicted proteins associated with all virophage markers detected, while
XXXX_virophage_cds_to_family.tsv includes the list of MCP hits used to assign new sequences to virophage families. All other intermediary
files can be safely ignored by most users.
