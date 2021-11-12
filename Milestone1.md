## Section 1: Update
* Download gdc-client for usage on command line following instructions for the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/#downloading-the-gdc-data-transfer-tool)
* Use gdc-client to download files from the [GDC Portal manifest files] (https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/)
	* Make sure manifest file is stored in an easily accessible directory 

```
gdc-client download -m  ~/final_project/gdc_manifest.specificfilenamehere.txt

```

* Use R studio to organize files into sigmoid and ascending colon cancer files (see in scripts)

![Example of Final Project Folder](/Screenshots/pic.png?raw=true)

![Example of Progress in DESeq2 Vignette](/Screenshots/pic2.png?raw=true)


## Section 2: Next Steps
* Complete initial draft running the entire vignette through htseq for feedback.
* Learn more about  HTSeq in python and DESeq2 packages and the differences
	* Check for any errors

## Section 3: Data
* Data will be stored in a GoogleDrive folder as a tar/zip file called ```finalproject.tar.gz.``` Usage of this data will ```tar -xvzf finalproject.tar.gz``` 
## Section 4: Known Issues 
* Retreiving data from GDC Portal through GCD Data Transfer Tool. ```gcd-client``` did not initially want to download using ```conda install```. I was able to work around this by installing via ubuntu and copying it to my $PATH. so I was unable to use it in the command line. 
