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