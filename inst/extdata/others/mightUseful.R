# prepare the txdb file
gr <- import("hg19_toy.gtf",asRangedData = FALSE)
hg19_toy <- makeTxDbFromGRanges(gr)
saveDb(hg19_toy, "hg19_toy.sqlite")
hg19_toy <- loadDb("hg19_toy.sqlite")

gr <- import("mm10_chr1_part.gtf")
mm10_toy <- makeTxDbFromGRanges(gr)
saveDb(mm10_toy, "mm10_toy.sqlite")
mm10_toy <- loadDb("mm10_toy.sqlite")

# package
package.skeleton(name="Travis",code_files =list.files())

# prepare other data sets
a <- read.table("cellular_dox_macs_kpdup_peaks.narrowPeak",sep="\t",header=FALSE)
write.table(a[1:1000,],file= "m6A_hg19_1000peaks_macs2.narrowPeak",sep="\t",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)

a <- read.table("mmDNA.bismark.cov",sep="\t",header=FALSE)
write.table(a[1:10000,],file= "DNAm5C_mm10_10000sites_bismark.cov",sep="\t",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)


a <- read.table("mmRNA.bismark.cov",sep="\t",header=FALSE)
write.table(a[1:10000,],file= "RNAm5C_mm10_10000sites_bismark.cov",sep="\t",
            quote = FALSE, row.names = FALSE,
            col.names = FALSE)


# prepare the bam file
samtools index 866991.bam
samtools index 866992.bam
samtools view -b 866991.bam chr1:6,177,566-7,207,283 -o 866991_part.bam
samtools view -b 866992.bam chr1:6,177,566-7,207,283 -o 866992_part.bam

7,236,448