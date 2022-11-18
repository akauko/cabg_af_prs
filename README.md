# cabg_af_prs
R code for the publication: Polygenic Risk Scores for Predicting Adverse Outcomes After Coronary Revascularization
https://doi.org/10.1016/j.amjcard.2021.11.046

* Data: FinnGen (https://www.finngen.fi/en)
* PRS values were calculuted for FinnGen individuals using PRS-CS pipeline with default settings: https://github.com/getian107/PRScs
* We used GWAS summaries from UKBB GWAS v3: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291

```
cabg_af_prs
├── README.md                 		# Overview
├── cabg_prs_after_download.rmd     	# R markdown for the analysis with all outcomes
├── cabg_af_prs_download.rmd     	# R markdown for the detailed analysis with AF only
├── functions2.R      			# Minor R functions for the main analysis
├── select_columns.pl         		# Perl script to select columns from tsv files by column name
├── prs_calculations		  		# Directory: PRS calculations
	├── preprocess_ukb_for_prscs.bash	# Preprocesses ukb gwas summary statistics for PRS-CS
	├── prs_ukb3.wdl			# Runs PRS-CS at FinnGen for one endpoint and all sexes.
	├── generate_prs_scores.bash		# Processes PRS_CS output and calculates weights for FinnGen individuals
	├── reformat_raw_prs.R 			# Is required by 'generate_prs_scores.bash'

```
