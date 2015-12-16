source("~/scripts/tags.r")
source("~/scripts/peaks.r")

#####

if(!interactive()) {
  fname <- commandArgs(T)[1]
}

peaks <- read.narrow.peak.file("peaks_peaks.narrowPeak")

pets <- readInBamFile(fname)
l <- length(pets)


joint <- !is.na(pets@pos) & !is.na(pets@mpos)

idx <- seq(1,sum(joint),by=2)


left <- GRanges(seqnames=as.character(pets@rname)[joint][idx],
              ranges=IRanges(pets@pos[joint][idx],width=1),
              strand='*')
right <- GRanges(seqnames=as.character(pets@rname)[joint][idx+1],
                ranges=IRanges(pets@pos[joint][idx+1],width=1),
                strand='*')

####

