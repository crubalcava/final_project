## Section 1: Milestone 1 Updates
* Download gdc-client for usage on command line following instructions for the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/#downloading-the-gdc-data-transfer-tool)
* Use gdc-client to download files from the [GDC Portal manifest files] (https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/)
	* Make sure manifest file is stored in an easily accessible directory 

```
gdc-client download -m  ~/final_project/gdc_manifest.specificfilenamehere.txt

```

* Use R studio to organize files into sigmoid and ascending colon cancer files (see in scripts)

![Example of Final Project Folder](/Screenshots/pic.png?raw=true)

![Example of Progress in DESeq2 Vignette](/Screenshots/pic2.png?raw=true)


## Section 2: Milestone 2 Updates
* Initial draft of DESeq2 vignette updated
* Use python (jupyterlab notebook) and R (RStudio) to fine tune file organization
	* htseq files successfully extracted from folders and organized into more effective sigmoid and ascending colon files
![htseq-counts.csv](/Screenshots/count_csv_format.png?raw=true)
	* manual workarounds using combine values to establish lists for sample tables in DESeq2 vignette
![Environment of DESeq2 Vignette](/Screenshots/milestone2_env.png?raw=true)

* Explore different result plots for analysis
	* [Volcano plot vignette](https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html#quick-start) 

## Section 3: Final Steps (DUE DECEMBER 3)
* Prepped count tables to join ascending and sigmoid specific ```htseq.counts``` files.
* Fix error messages found in [DESeq2Update](https://github.com/crubalcava/final_project/blob/main/Scripts/DESeq2Update.Rmd) 
	* Error message from following the [DESeq2 vignette] (http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq-count-input)
* Run graphs and analyze for differential gene expression.
* Clean up vignette and combine final [scripts](https://github.com/crubalcava/final_project/tree/main/Scripts) to properly guide others.
* change ens id to gene id to make data more comprehensive
	
## Section 4: Data
* Data will be stored in a GoogleDrive folder as a tar/zip file called ```finalproject.tar.gz.``` Usage of this data will ```tar -xvzf finalproject.tar.gz``` 
## Section 5: Known Issues 
* Retreiving data from GDC Portal through GCD Data Transfer Tool. ```gcd-client``` did not initially want to download using ```conda install```. I was able to work around this by installing via ubuntu and copying it to my $PATH. so I was unable to use it in the command line. 
* Multiple error messages of ```Error in file(file, "rt") : cannot open the connection``` or ```Warning in file(file, "rt") :  cannot open file...No such file or directory ```
	* avoid this error by setting working directory. may need to be done multiple times in R Notebook OR you can assign a working directory when creating the notebook. 