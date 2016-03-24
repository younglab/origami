---
output: pdf_document
---
# origami

## Installation

### Depedencies
* cutadapt (version 1.8)
* samtools (version 1.2)
* bamtools
* bowtie
* MACS2
* R

## Running the pipeline

The pipeline has two core steps that need to be run:

* origami-alignment -- aligns the ChIA-PET reads
* origami-analysis -- estimates the confidence in each ChIA-PET interaction

### Aligning ChIA-PET reads (origami-alignment)

Aligning the reads is handled through the origami-alignment executable.

Command: origami-alignment [options] <bowtie idx> <FASTQ read 1> <FASTQ read 2>

This script will take a bowtie index (*must* be Bowtie 1 index) and the paired read files and handle the adapter trimming, alignment, and read peak calling. The output of the executable will be a BAM file of aligned reads as well as a list of peaks called (via MACS2). The FASTQ files can be in either gzipped or bz2zipped formats and the software will correctly unpack them. The optional arguments are below:

* -o,--output=: The output directory to put all the aligned reads and peak calls in. Defaults to "output"
* -h: Shows the help menu
* -v: Turn on verbose mode
* -p: Activates parallization on the LSF queue, *must* be paired with --lsf-queue. It is also recommended that --splitnum be set approrpiate to the data
* -m=,--min-len=[integer]: Sets the minimum read length to keep post-adapter trimming. Default is 15 bp
* --keep-tmp: Don't delete temporary files when the alignment is finished
* --lsf-queue=[queue name]: Sets the LSF queue to queue all the parallized jobs, implies -p. There is no default for this option
* --forward-linker: Set the forward linker to search for in the paired reads. Default is ACGCGATATCTTATCTGACT.
* --reverse-linker: Set the reverse linker to search for in the paired reads. Default is AGTCAGATAAGATATCGCGT.
* --ab-linker: Activates the AB-linker trimming mode instead of the default long-read linker mode.
* --a-linker: Sets the A-linker sequence. Default is CTGCTGTCCG. Implies --ab-linker
* --b-linker: Sets the B-linker sequence. Default is CTGCTGTCAT. Implies --ab-linker
* --pp=[executable command]: Preprocess reads. Runs each of the passed FASTQ files through the executable script before trimming and alignment. Default is none
* --macs-gsize=[value]: Set the genome size option for MACS2 (see MACS2 for how this option is set). Default is hs


### Estimating the confidence estimates in each interaction (origami-analysis)



Command: origami-analysis [options] <BAM file> <peaks file> <output prefix>

After the alignment step, the next executable estimates the confidence in each identified interaction within the data set. This takes the resultsing BAM file and peak file from the alignment step and handles the identification of interactions and estimating their confidence. Any BAM file that is paired end and shorted by read name can be passed to the script. Any BED-type peaks file can be passed to the algorithm. The executable will generate two files: <output prefix>-model-data.Rdata and <output prefix>results.csv. The options are below:

* -h: Shows the help menu
* --without-distance: Turns off weighting interactions by their ditsance
* --iterations=[positive integer]: Sets the number of iterations for the MCMC algorithm. Default is 10000
* --burn-in=[integer]: Sets the number of iterations used for the burn-in step of the MCMC algorithm. Default is 100. A value of 0 turns off the burn-in period.
* --prune=[0 or 2+]: Sets the number of steps at which the MCMC algorithm prunes an interation. The default is 5. A value of 0 turns off the pruning.
* --slop-dist=[integer]: Extends peak peak in both directions by the given length. Default is 0. This also requires --slop-genome
* --slop-genome=[file]: A file that contains two columns: the first for each chromosome in the genome and the length of the chromosome. Needs --slop-dist.
* --save-full-model: Saves all the MCMC data generated. By default it saves a reduced set of essential parameters. (*Warning*: these files can get very big and can potentially exhaust the RAM while running, so use this for a low number of iterations if possible)

### Converting the results to a useable format (origami-conversion)


From the results CSV file generated from the origami-analysis, the results file can be converted through the following command:

Command: origami-conversion <bedpe/bed/washu> <results CSV file> <column index>

This will convert the interactions in the results file with their particular column into a BEDPE, BED12, or WashU file format. For the BED12 file format, since this assumes each entry starts and ends on the same chromosome, this excludes all interchromosmal interations. Also, this will report all interactions regardless of significance, so one may want to filter the results after the conversion. The script allows for the following column indexes (generally one will use 3):

* 1 -- Reports the PET count of each interaction
* 2 -- Reportes the (uncorrected) hypergeometric p-value of the interaction (based on the Fullwood et al 2009 method)
* 3 -- Reports the confidence/belief probability from the MCMC procedue