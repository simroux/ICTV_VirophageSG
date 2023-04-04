# ICTV_VirophageSG
Data and information related to work by the ICTV Virophage Study Group

# Databases

- All_markers.hmm  
HMM database of virophage-specific markers intended to assign new genomes to the *Maveriviricetes* class.

- MCP_blast_db  
Blast database intended to assign members of the *Maveriviricetes* class to lower taxa. Note: ideally, this assignment would be based on a phylogenomic analysis.
This blast database is only provided to offer an alternative method to assign new sequences to existing virophage taxa based, in cases where a de novo phylogenomic
analysis is impractical.

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

## Genus delineation and example scripts

Genera within the *Maveriviricetes* class are proposed to be defined based on the identification of 7 shared genes or more between pairs of genomes. An example of this gene-sharing-based clustering is provided in "Example_scripts_genus/"

# About - Copyright Notice

Custom script used in 'Updated virophage taxonomy and distinction from
polinton-like viruses' (Virophage_Scripts) Copyright (c) 2023, The Regents
of the University of California, through Lawrence Berkeley National Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of Energy). 
All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Intellectual Property Office at
IPO@lbl.gov.

NOTICE.  This Software was developed under funding from the U.S. Department
of Energy and the U.S. Government consequently retains certain rights.  As
such, the U.S. Government has been granted for itself and others acting on
its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the
Software to reproduce, distribute copies to the public, prepare derivative 
works, and perform publicly and display publicly, and to permit others to do so.


