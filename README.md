# SIV-R01

## Analysis of SIV R01 SIV/Mtb Co-infection Study


* Data directory contains sample meta data/information. SampleSheet_Mac251.csv has the most comprehensive sample information.
* Output directory contains all output files, except the bam files (due to the memory limit of GitHub).
* Rscripts directory contains all scripts for the analysis. 


## Analysis steps

### FASTQ files are processed by the PID (Primer ID) pipeline

### Step 1. Rename the fastq directory to two digit number (1_ or 10 etc.) 
* Run the first part of '1.CreateBash_PID.R'
* This will create a bash script 'rename.RunX.sh' in the directory where FASTQ are located

```bash
#Run it by
bash rename.RunX.sh
```

### Step 2. Create bash files to run the PID pipeline from raw FASTQ files
* Use the 'Step 2' section of 1.CreateBash_PID.R
* This will create bash scripts under PID-master directory '~/programs/PID-master/Bashscripts/'
* Run all scripts in PID-master directory 
  - The bash file will run BBTools to trim low quality bases from a pair of raw fastq files and run the PID pipeline on them.
  - The output files from PID are mapped to the reference sequence using BWA.
  - Sam files are hard-clipped for merging forward and reverse reads in R (Step 3)

```bash
#Run all bash scripts in a directory
for f in BashScripts/*.sh; do
	bash "$f" -H   || break 
done
```
* OUTPUT: Hard-clipped sam files in 'PID-master/Sam/RUNX_XX.clipped.sam'


### Step 3. Merge forward and reverse reads from Sam files using MergeForRevPID.R. 
* Output is in PID_Consensus/ as a fasta file (one per sample).

### Step 4. Map merged sequences (from Step 3) to the reference genome using bwa with a relaxed setting 
* Use the 'Step 3' section of 1.CreateBash_PID.R to create all bash files
* Run all of them. Output is sorted bam files and bam index files, saved in Output/bam_PID/

	
## From here, you can conduct either diversity/divergence (mutation frequency) analysis (Step 5) or focus on interesting (high-frequency) sites and analyze 'genotype' data 
  
### Step 5. Use the pileup to create mutation frequency tables and conduct analyses. Use... 
* Scripts 2 to 4.1 create frequency tables from each bam file.
* Scripts 4.1.2 is used to remove the center portion from Run2 files
* Script 4.2 to assess indels 
* Script 4.3 to look at the stock and control files
* Script 5.1 to look at mutation frequency/diversity patterns of each sample and within each monkey
* Script 5.2 to find sites with a high mutation frequency to select the sites of interest for further analyses     
* Script 5.3 compares the mutation frequencies of samples with and without Primer ID
* Script 5.4 compares viral load and viral diversity 
    
    
### Step 6. Create genotype data/figures using 6.0/6.1 scripts for nucleotide, and 7.0-7.5 for amino acids
* Script 6.0 to create an alignend fasta file of consensus PID sequences for each file
* Script 6.1 to create 'genotype' data based on nucleotides using the sites identified by Script 5.2 (currently 11 sites)
* Script 7.0 to create fasta files with translated amino acids for each sample 
* Script 7.1 to create 'genotype' data based on amino acids
* Script 7.2 to create genotype figures (stacked area charts, streamgraphs)
* Script 7.3 to assess how AA changed over time in each monkey  
* Script 7.4 to focus on the 3 mutations (D99E, K102R, P144Q) -create genotype freq based on the 3 sties
* Script 7.5 to further look into the 3 mutations  


### Other Analyses
* Script Genotype_linkage.R: Explore linkage between sites  
* Script Haplotype_plots.R: Create figures of aligned (haplotype) sequences (50-60 sequences/sample)
* Script subSample_reads.R: Create figures of variable sites from sub-sampled reads (similar to Haplotype_plots, but more focused, old analysis)
* Script Compare_PID_withTrimmedPID.R: Raw reads were trimmed at Q15 before processing by PID pipeline.
* Script MF_PIDcon_compare.R: Compare PID processed files without creating consensus sequences vs. with consensus reads (merging Forward and Reverse to create each consensus read)
* Script LD_plots.R: Create figures based on LD data from CoVaMa 
	
###  Scripts starting with A are for samples without PID
