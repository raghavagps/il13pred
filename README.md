# **IL13Pred**
## Introduction
IL13Pred is developed for predicting, desiging, and scanning the interleukin-13 inducing peptides. More information on IL13Pred is available from its web-server https://webs.iiitd.edu.in/raghava/il13pred/stand.html. This page provides information about stnadalone version of IL13Pred. Please read/cite the content about the IL13Pred for complete information including algorithm behind IL13Pred.

* Models: In this program, one model has been incorporated for predicting interleukin-13 inducing peptides. The model is trained on IL-13 inducing and non-inducing peptides.

* Modules/Jobs: This program implements three modules (job types); i) Predict: for predictin of interleukin-13 inducing peptides, ii) Design: for generating all possible mutant peptides and computing interleukin-13 inducing potential (score) of peptides, iii) Scan: for creating all possible overlapping peptides of given length (window) and computing interleukin-13 inducing potential (score) of these overlapping peptides.

* Minimum USAGE: Minimum usage is "python il13pred.py -i peptide.fa," where peptide.fa is a input fasta file. This will predict the interleukin-13 inducing potential of sequence  in fasta format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

* Full Usage: Following is complete list of all options, you may get these options by "python il13pred.py -h" 
