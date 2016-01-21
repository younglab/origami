source("~/scripts/tags.r")
source("~/scripts/peaks.r")

read.bed.file <- function(file) {
  temp <- read.table(file,sep='\t')
  
  GRanges(seqnames=as.character(temp$V1),ranges=IRanges(temp$V2,temp$V3),strand='*')
}

if(!interactive()) {
  args <- commandArgs(T)
  bamfile <- if(is.na(args[1])) "mapped_reads.bam" else args[1]
  peakfile <- if(is.na(args[2])) "peaks_peaks.narrowPeak" else args[2]
  peakcountsfile <- if(is.na(args[3])) "peak-counts.txt" else args[3]
  intcountsfile <- if(is.na(args[4])) "int-counts.txt" else args[4]
  
}

if( is.null(bamfile) || is.null(peakfile) ) {
	stop("estimate-interaction-counts.r <BAM file> <peaks file> [peak counts file] [interaction counts file]")
}

pets <- readInBamFile(bamfile)
peaks <- read.bed.file(peakfile)

both.mapped <- !is.na(pets@pos) & !is.na(pets@mpos)

idx <- seq(1,sum(both.mapped),by=2)

first <- GRanges(seqnames=pets@rname[both.mapped][idx],
                 ranges=IRanges(pets@pos[both.mapped][idx],width=1),
                 strand='*')
second <- GRanges(seqnames=pets@rname[both.mapped][idx+1],
                  ranges=IRanges(pets@pos[both.mapped][idx+1],width=1),
                  strand='*')

of <- findOverlaps(first,peaks)
os <- findOverlaps(second,peaks)

total <- table(c(subjectHits(of),subjectHits(os)))
p <- peaks
mcols(p)[,"counts"] <- rep(0,length(p))
mcols(p)[,"counts"][as.integer(names(total))] <- as.vector(total)

write.table(as.data.frame(p)[,c(1,2,3,4)],file=peakcountsfile,sep='\t',col.names = F,row.names = F,quote=F)

l <- list()
qof <- queryHits(of)
sof <- subjectHits(of)
qos <- queryHits(os)
sos <- subjectHits(os)

midx <- match(qof,qos)
m <- matrix(c(sof[!is.na(midx)],sos[midx[!is.na(midx)]]),ncol=2)
m <- m[order(m[,1],m[,2]),]

f <- factor(paste(m[,1],m[,2],sep='_'))

intcounts <- do.call(rbind,lapply(split(1:nrow(m),f),function(idx) {
  if(length(idx) == 1 ) {
    ret <- c(m[idx,],1)
  } else {
    ret <- c(m[idx[1],],length(idx))
  }
  
  ret
}))

rownames(intcounts) <- NULL

outtable <- as.data.frame(peaks)[,c(1,2,3)]

locations <- cbind(outtable[intcounts[,1],],outtable[intcounts[,2],])

samechromosome <- locations[,1] == locations[,4]

if(any(samechromosome)) {
  flip <- as.integer(outtable[intcounts[,1],][samechromosome,2]) > as.integer(outtable[intcounts[,2],][samechromosome,2])
  outtable <- cbind(outtable[intcounts[,1],],outtable[intcounts[,2],],intcounts[,3])
  
  
  tmp <- outtable[samechromosome,][flip,4:6]
  outtable[samechromosome,][flip,4:6] <- outtable[samechromosome,][flip,1:3]
  outtable[samechromosome,][flip,1:3] <- tmp
} else {
  outtable <- cbind(outtable[intcounts[,1],],outtable[intcounts[,2],],intcounts[,3])
}
write.table(outtable,file=intcountsfile,sep='\t',col.names=F,row.names=F,quote=F)
