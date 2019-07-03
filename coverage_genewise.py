import os
import subprocess as sp
import pandas as pd

# indir = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Bams/"
# wd = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/"
# os.chdir(wd)

# bams = [b for b in os.listdir(indir) if b.endswith("q5.bam")]

# # 1. Remove Duplicates

# for bam in bams:
#     i = indir+bam
#     o = wd+bam.replace(".bam", "_noDup.bam")
#     m = wd+bam.replace(".bam", "_metrics.txt")

#     cmd = ("java -jar $PICARD MarkDuplicates "
#            "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

#     sp.call(cmd, shell=True)

# # 2. Index new bams

# bams = [b for b in os.listdir(wd) if b.endswith("_noDup.bam")]

# for bam in bams:
#     cmd = "samtools index {}" .format(bam)
#     sp.call(cmd, shell=True)

# # 3. Calculate number of reads

# nreads = {}
# for b in bams:

#     cmd = "samtools view -c {}" .format(b)
#     nreads[b] = int(sp.check_output(cmd, shell=True))

# # 4. Calculate coverage


# def coverage(bam, ref, out):
#     cmd = "bedtools multicov -bams {} -bed {} -p > {}" .format(bam, ref, out)
#     sp.call(cmd, shell=True)


# bamlist = " ".join(bams)
# ref = "/home/lucas/ISGlobal/Gen_Referencies/PlasmoDB-43-gene_ref.bed"
# coverage(bamlist, ref, "all_coverage.bed")

# ref_prom = "/home/lucas/ISGlobal/Gen_Referencies/promoters2000bp_ref.bed"
# coverage(bamlist, ref_prom, "all_coveragePromoters.bed")


# # 5. Normalize by Read Count and kb

# def normalizeRPKM(bed):
#     header = ["chrom", "start", "stop", "gene_id"] + bams

#     with open(bed, "r+") as infile:
#         with open(bed.replace(".bed", "_RPKM.bed"), "w+") as outfile:
#             outfile.write("\t".join(header)+"\n")
#             for line in infile:
#                 linelist = line.strip().split("\t")
#                 raw_vals = [int(v) for v in linelist[4:]]
#                 genlen = int(linelist[2]) - int(linelist[1])
#                 rpkm_vals = []
#                 for i in range(len(bams)):
#                     rpm = raw_vals[i]/(nreads[bams[i]]/2000000)
#                     rpkm = rpm/(genlen/1000)
#                     rpkm_vals.append(rpkm)

#                 rpkm_str = [str(v) for v in rpkm_vals]

#                 outfile.write("\t".join(linelist[0:4]+rpkm_str)+"\n")


# normalizeRPKM("all_coverage.bed")
# normalizeRPKM("all_coveragePromoters.bed")

# 6. Normalize by Input

wd = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/GenWise_Coverage/"
os.chdir(wd)

rpkm_beds = [b for b in os.listdir(wd) if b.endswith("_RPKM.bed")]
rpkm_beds

# Substract Input

df = pd.read_csv(rpkm_beds[2], sep="\t")
df

for col in df.columns[4:]:
    mark = col.split("_")[1]
    prfx = col.split("_")[0]
    if mark in ["me", "ac"]:
        treat = col
        ctl = [x for x in df.columns[4:]
               if x.split("_")[0] == prfx and
               x.split("_")[1] == "in"][0]

        df[prfx] = df[treat] - df[ctl]

# df.to_csv("prova_normIN.csv")
