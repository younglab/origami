source("~/scripts/tags.r")
source("~/scripts/peaks.r")

if(!interactive()) {
  args <- commandArgs(T)
  bamfile <- args[1]
  peakfile <- args[2]
}

pets <- readInBamFile(bamfile)
peaks <- read.narrow.peak.file(peakfile)

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

write.table(as.data.frame(p)[,c(1,2,3,12)],file='peak-counts.txt',sep='\t',col.names = F,row.names = F,quote=F)

l <- list()

for( i in 1:length(p) ) {
  l[[i]] <- unique(c(queryHits(of)[subjectHits(of)==i],queryHits(os)[subjectHits(os)==i]))

}
