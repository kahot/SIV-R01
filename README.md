# SIV-R01

## Analysis of SIV R01 SIV/Mtb Co-infection Study


* Data directory contains sample meta data/information. SampleSheet_Mac251.csv has the most comprehensive sample information.
* Output directory contains all output files, except the bam files (due to the memory limit of GitHub).
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

	
# From here, you can conduct either diversity/divergence (mutation frequency) analysis (Step 5) 
#  or focus on interesting (high-frequency) sites and analyze 'genotype' data 
  
### Step 5. Use the pileup to create mutation frequency tables and conduct analyses. Use... 
    • Scripts 2 to 4.1 create frequency tables from each bam file.
    • Script 4.2 to assess indels 
    • Script 4.3 to look at the stock and control files
    • Script 5.1 to look at mutation frequency/diversity patterns of each sample and within each monkey
    • Script 5.2 to find sites with a high mutation freqeuncy to select the sites of interest for further analyses     
    • Script 5.3 compares the mutation frequenceis of samples with and without Primer ID
    
    
### Step 6.1 Create genotype tables using 6+ scripts.
	• Script 6.0 to create an aligend fasta file of consensus PID sequences for each file
	• Script 6.1 to create 'genotype' data bassed on nucleotides using the sites identified by Script 5.2 (currently 11 sites)
	• Script 7.0 to create fasta files with translated amino acids for each sample 
	• Script 7.1 to create 'genotype' data bassed on amino acids
	• Script 7.2 to create genotype figures (stacked area charts, streamgraphs)
	• Script 7.3 to assess how AA changed over time in each monkey  
    • Script 7.4 to focus on the 3 mutations (D99E, K102R, P144Q) -create genotype freq based on the 3 sties
    • Script 7.5 to further look into the 3 mutations  
    

###  Scripts starting with A are for samples without PID
