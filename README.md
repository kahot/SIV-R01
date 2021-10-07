# SIV-R01

## Analysis of SIV R01 SIV/Mtb Co-infection Study


* Data directory contains sample information. SampleSheet_Mac251.csv has the most comprehensive sample information.
* Output directory contains all output files, except the bam files (due to the memory limit)
* Rscripts directory contains all scripts for the analysis. 


## Analysis steps

### FASTQ files are processed by the PID (Primer ID) pipeline

### Step 1. Rename the fastq directory to two digit number (1_ or 10 etc.) 
	• Run the first part of '1.CreateBash_PID.R'
	• This will create a bash script 'rename.RunX.sh' in the directory where FASTQ are located

```bash
#Run it by
bash rename.RunX.sh
```

### Step 2. Create bash files to run the PID pipeline from raw FASTQ files
	• Use the 'Step 2' section of 1.CreateBash_PID.R
	• This will create bash scripts under PID-master directory '~/programs/PID-master/Bashscripts/'
	• Run all scripts in PID-master directory

```bash
#Run all bash scripts in a directory
for f in BashScripts/*.sh; do
	bash "$f" -H   || break 
done
```
	• Hard-clipped sam files will be saved in 'PID-master/Sam/RUNX_XX.clipped.sam'


### Step 3. Merge forward and reverse reads from Sam files using MergeForRevPID.R. 
	• Output is in PID_Consensus/ as a fasta file (one per sample).

### Step 4. Map merged sequences (from Step 3)to the reference genome using bwa with a relaxed setting 
	• Use the 'Step 3' section of 1.CreateBash_PID.R to create all bash files
	• Run all of them. Output is sorted bam files and bam index files, saved in Output/bam_PID/
	
### Step 5. A) Use the pileup to create A) frequency tables (2.pileup_Rsamtools_PID.R) or
###         B) Create  genotype data using 8.1.Genotyping.R  Output is in Output/Genotype/
	

### Analysis 5A: calculate mutation frequencies and diversity from frequency tables
    • Follow the R scripts #2-7 to assess diversity and divergence.  
	
	
#4. Create bash files to map PID-consensus.fasta to the reference using bwa with a relaxed setting. The template is Bash_mapConPID.sh. Run all bash files    
### 2. After mapping is done, start the analysis in R by following the file numbers

    • Scripts 1 to 3 create frequency tables for each bam file.
    • Script 4 analyze the stock and control files 
    • Script 5 assess trasnmitted/founder variants 
    • Script 6 assess SIV diversity across tissues/time/cohorts and compare with RNA copy nubmers and CD4/CD8 frequency
    • Script 7 analyze the indels 
    • Script 8 assess the known immune escape mutations 
    • Script 9 assess an intereting mutation at AA112
    • Script 10 analyze the effects of genetic drift between tissues
    • Script 11 conducts glmm to udnerstand the effects of different factors on SIV diversity
    
