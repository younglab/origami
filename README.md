# origami

## Overview

Origami is a pipeline for processing and calling high-confidence chromatin loops associated with 
the ChIPped factor. The pipeline implements a semi-Bayesian two-compenent mixture model to learn
which contacts within the data set appear to represent apparent structured chromatin loops
in a cell population while correcting for important biases in ChIA-PET data (notably, linear
genomic distance between the ends of a potential chromatin loop).

## Contact Information

For more information or any questions, please contact dsday@wi.mit.edu or young_computation@wi.mit.edu.

## Installation

To install origami, first run the configure script to test for the presence of mdost of the dependencies origami needs. After configure completes successfully, run make to compile the C++ code. Then  make install will copy all the necessary files to the installation directory.

(Note: at the present time, the configure script does not check if Bamtools is installed correctly, please make sure it is! The Makfile may need to be adjusted to adapt to where Bamtools API files are installed on your system.)

### Depedencies
* cutadapt (tested with version 1.8)
* samtools (tested with version 1.2)
* bamtools
* bowtie1
* MACS1 and MACS2
* R (also requires the R libraries GenomicsRanges, matrixStats)
* g++
* bash
* make
* autoconf


## Running the pipeline

The pipeline has three steps that need to be run:

* origami-alignment -- aligns the ChIA-PET reads
* origami-analysis -- estimates the confidence in each ChIA-PET interaction
* origami-conversion -- converts the origami output to popular genomic file formats for visualization

### Aligning ChIA-PET reads (origami-alignment)

Aligning the reads is handled through the origami-alignment executable.

Command: origami-alignment [options] <bowtie idx> <FASTQ read 1> <FASTQ read 2>

This script will take a bowtie index (*must* be Bowtie 1 index) and the paired read files and handle the adapter trimming, alignment, and read peak calling. The output of this step is a BAM file of aligned paired reads and peak calls from MACS1 and MACS2 named accordingly. The FASTQ files can be in either gzipped or bz2zipped formats. The optional arguments are below:

* -o,--output=: The output directory to put all the aligned reads and peak calls in. Defaults to "output"
* -h: Shows the help menu
* -v: Turn on verbose mode
* -p: Activates parallization on the LSF queue, *must* be paired with --lsf-queue. It is also recommended that --splitnum be set approrpiate to the data
* --divide-pets=[NUM]: When -p is active, this will divide the FASTQ reads into subfiles of NUM reads to distribute over a LSF queue for faster processing
* --lsf-queue=[queue name]: Sets the LSF queue to queue all the parallized jobs, implies -p. There is no default for this option
* -m=,--min-len=[integer]: Sets the minimum read length to keep post-adapter trimming. Default is 15 bp
* --mode=[mode]: Sets the type of linker-trimming mode, either long (bridge-linker) or ab (AB-linker) are acceptable options
* --forward-linker: Set the forward linker to search for in the paired reads. Default is ACGCGATATCTTATCTGACT.
* --reverse-linker: Set the reverse linker to search for in the paired reads. Default is AGTCAGATAAGATATCGCGT.
* --ab-linker: Activates the AB-linker trimming mode instead of the default long-read linker mode.
* --a-linker: Sets the A-linker sequence. Default is CTGCTGTCCG. Implies --ab-linker
* --b-linker: Sets the B-linker sequence. Default is CTGCTGTCAT. Implies --ab-linker
* --pp=[executable command]: Preprocess reads. Runs each of the passed FASTQ files through the executable script before trimming and alignment. Default is no preprocessing
* --macs-gsize=[value]: Set the genome size option for MACS1 and MACS2 (see the MACS documentation). Default is hs
* --keep-tmp: Don't delete temporary files when the alignment is finished



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
