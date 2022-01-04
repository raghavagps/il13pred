# **IL13Pred**
## Introduction
IL13Pred is developed for predicting, desiging, and scanning the interleukin-13 inducing peptides. More information on IL13Pred is available from its web-server https://webs.iiitd.edu.in/raghava/il13pred/stand.html. This page provides information about stnadalone version of IL13Pred. Please read/cite the content about the IL13Pred for complete information including algorithm behind IL13Pred.

**Models:** In this program, one model has been incorporated for predicting interleukin-13 inducing peptides. The model is trained on IL-13 inducing and non-inducing peptides.

**Modules/Jobs:** This program implements three modules (job types); i) Predict: for predictin of interleukin-13 inducing peptides, ii) Design: for generating all possible mutant peptides and computing interleukin-13 inducing potential (score) of peptides, iii) Scan: for creating all possible overlapping peptides of given length (window) and computing interleukin-13 inducing potential (score) of these overlapping peptides.

**Minimum USAGE:** Minimum usage is "python il13pred.py -i peptide.fa," where peptide.fa is a input fasta file. This will predict the interleukin-13 inducing potential of sequence  in fasta format. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

**Full Usage:** Following is complete list of all options, you may get these options by "python il13pred.py -h" 

***Usage:*** *il13pred.py [-h] -i INPUT* 
		<br>    *[-o OUTPUT]*
		<br>	*[-j {1,2,3}]*
		<br>	*[-t THRESHOLD]* 
		<br>	*[-w {8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35}]* 
		<br>	*[-d {1,2}]*

***Arguments Description:***
 <br> -h,   *--help: show help message and exit.*
 <br> -i INPUT, *--input: protein or peptide sequence in FASTA format or single sequence per line in single letter code.*
  <br> -o OUTPUT, *--output: File for saving results by default outfile.csv.*
  <br> -j, *--Job Type: 1:predict, 2:design and 3:scan, by default 1.*
 <br> -t THRESHOLD, *--Threshold: Value between 0 to 1 by default 0.06.*
 <br> -w, or --winleng, *--Window Length: 8 to 35 (scan mode only), by default 9.*
  <br> -d, *--Display: 1:Interleukin-13 inducing peptide, 2: All peptides, by default 1*


**Input File:** It allow users to provide input in two format; i) FASTA format (standard) and ii) Simple Format. In case of simple format, file should have one peptide sequence in a single line in single letter code (eg. peptide.seq). 


**Note:**
<br> 1: In case of predict and design module (job), the length of peptide should be upto 35 amino acids. If a sequence with length more than 35 will be provided, program willtake first 35 residues, and ignore the rest. In case of scan module, minimum length of protein/peptide sequence should be more than or equal to window length (pattern), see peptide.fa.
<br> 2: Program will ignore peptides having length less than 8 residues (e.g., protein.fa).

**Output File:** Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

**Threshold:** User should provide threshold between 0 and 1, please note that the score is propotional to interleukin-13 inducing potential of peptide.


## IL13Pred Package Files
=======================
<br> Brief description of the files included is given below:

* INSTALLATION  			: Installations instructions

* LICENSE       			: License information

* README.md     			: This file provide information about this package

* XGB_model       		: Model file comprising the parameters of XGB classifier

* il13pred.py 			: Main python program 

* peptide.fa			: Example file contain peptide sequenaces in FASTA format

* peptide.seq			: Example file contain peptide sequenaces in simple format

* protein.fa			: Example file contain protein sequenaces in FASTA format 

* example_predict_output.csv	: Example output file for predict module

* example_scan_output.csv		: Example output file for scan module

* example_design_output.csv	: Example output file for design module
            	
