# Final Project Outline
## Title 
Differential Gene Expression in TCGA within Stage II Colon Adenocarcinomas comparing ascending colon and sigmoid colon tumor orgin sites using DeSEQ2

## Author 
Crystal Rubalcava

## Overview of Project 
I will identify differentially expressed genes between Colon Cancer Adenocarcinomas for patients with tumor sites originating in the ascending colon versus the sigmoid colon. Tumors found in the ascending colon are usually associated with right-sided colorectal (colon) cancer. Tumors found in the sigmoid colon are typically associated with left-sided colorectal (colon) cancer. This is important because  colorectal (colon) cancer does not contain a single type of tumor and origin site of the tumor can determine the molecular characteristics of the disease. To learn more about this difference, you can refer to the literature review here: [Difference Between Left-Sided and Right-Sided Colorectal Cancer: A Focused Review of Literature.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6089587/)

This analysis will utilize the package DESeq2 and follow the specific vignette below for htseq-count files.

* Vignette:[http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html.](http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 

For this analysis, I'll use the TCGA cohort and have identified 85 ht-seq counts files for tumors that fit within my cohort with 36 cases of patients with tumor sites originating in the ascending colon and 33 cases of patients with tumor sites originating in the sigmoid colon. I will control for the stage of the cancer and the disease type. 

## Data
I will use the clinical data from [portal.gdc.cancer.gov/repository](https://portal.gdc.cancer.gov/repository.) From the data available for colon adenocarcinomas are in ```htseq.count``` files. There are 85 tumor samples where 36 had tumor origin sites in the ascending colon and 33 had cases in sigmoid colon.


## Milestone 1
* Download fully loaded TCGS clinical data into a vignette through HTSeq steps.
* Start an initial draft of the analysis needed for the vignette.


## Milestone 2
* Complete initial draft running the entire vignette through htseq for feedback.


## Deliverable 
* A complete repository with clear documentation and description of your analysis and results.


