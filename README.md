# ICTV_VirophageSG
Data and information related to work by the ICTV Virophage Study Group

# Maveriviricetes_hallmark.2022-05-19.hmm
HMM models for the 3 Maveriviricetes hallmark genes: Major Capsid protein (MCP), ATPase, and Cysteine protase (Cprot).
Example of search using these models:
$ hmmsearch --cut_tc --noali -o Your_contigs-vs-hallmark.tc.full_output --tblout Your_contigs-vs-hallmark.tc.tsv Maveriviricetes_hallmark.2022-05-19.hmm  Your_contigs.fastaa
