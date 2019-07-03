import os
import subprocess as sp
import pandas as pd

indir = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/Renamed/RPM/"
wd = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/Renamed/RPM/"
os.chdir(wd)

""" All files must be in the same directory. 
All files must follow the naming convention: name_in/me/ac_... """


bams = [b for b in os.listdir(indir) if b.endswith("q5.bam")]
bams = [bams[0]]
bams

# 1. Remove Duplicates

for bam in bams:
    i = indir+bam
    o = wd+bam.replace(".bam", "_noDup.bam")
    m = wd+bam.replace(".bam", "_metrics.txt")

    cmd = ("java -jar $PICARD MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)

2. Index new bams

bams = [b for b in os.listdir(wd) if b.endswith("_noDup.bam")]
bams
# for bam in bams:
#     cmd = "samtools index {}" .format(bam)
#     sp.call(cmd, shell=True)

# 3. Calculate number of reads

nreads = {}
for b in bams:

    cmd = "samtools view -c {}" .format(b)
    nreads[b] = int(sp.check_output(cmd, shell=True))

# 4. Calculate coverage


def coverage(bam, out):
    cmd = "bedtools genomecov -ibam {} -pc -d > {}" .format(bam, out)
    sp.call(cmd, shell=True)


for bam in bams:
    coverage(bam, bam.replace(".bam", "_depth.prebed"))

# 5. Reformat depth file to proper BED

prebeds = [b for b in os.listdir(wd) if b.endswith("_depth.prebed")]


def properBED(prebed):
    out = prebed.replace(".prebed", ".bdg")
    with open(out, "w+") as outfile:
        with open(prebed, "r+") as infile:
            for line in infile:
                linelist = line.strip().split()
                chrom = linelist[0]
                start = str(int(linelist[1])-1)
                stop = linelist[1]
                depth = linelist[2]
                bedline = [chrom, start, stop, depth]
                outfile.write("\t".join(bedline)+"\n")


for pre in prebeds:
    properBED(pre)

# 6. Normalize by Number of Reads

beds = [b for b in os.listdir(wd) if b.endswith("_depth.bdg")]


def normalizeRPM(bed):

    with open(bed, "r+") as infile:
        with open(bed.replace(".bdg", "_RPM.bdg"), "w+") as outfile:
            for line in infile:
                linelist = line.strip().split("\t")
                raw_val = int(linelist[3])
                genlen = int(linelist[2]) - int(linelist[1])

                bamname = bed.replace("_depth.bdg", ".bam")
                rpm = raw_val/(nreads[bamname]/2000000)

                outfile.write("\t".join(linelist[0:3])+"\t"+str(rpm)+"\n")


for bed in beds:
    normalizeRPM(bed)

# 7. Normalize by Input (RPM/RPM)

rpm_beds = [b for b in os.listdir(wd) if b.endswith("_RPM.bdg")]
rpm_beds
# Make tupples met-in
# Filenames must be in the form: name_in/me/ac_...

# Divide by input
pseudo = 0.1

for b in rpm_beds:
    mark = b.split("_")[1]
    prfx = b.split("_")[0]
    if mark == "me" or mark == "ac":
        treat_file = b
        ctl_file = [x for x in rpm_beds
                    if x.split("_")[0] == prfx and
                    x.split("_")[1] == "in"][0]

        treat = pd.read_csv(treat_file, sep="\t", header=None)
        ctl = pd.read_csv(ctl_file, sep="\t", header=None)

        treat.iloc[:, 3] = (treat.iloc[:, 3]+pseudo)/(ctl.iloc[:, 3]+pseudo)
        treat.to_csv(path_or_buf=b.replace(
            ".bdg", "_normIN.bdg"), sep="\t", index=False, header=None)

# Substract Input
for b in rpm_beds:
    mark = b.split("_")[1]
    prfx = b.split("_")[0]
    if mark == "me" or mark == "ac":
        treat_file = b
        ctl_file = [x for x in rpm_beds
                    if x.split("_")[0] == prfx and
                    x.split("_")[1] == "in"][0]

        treat = pd.read_csv(treat_file, sep="\t", header=None)
        ctl = pd.read_csv(ctl_file, sep="\t", header=None)

        treat.iloc[:, 3] = treat.iloc[:, 3]-ctl.iloc[:, 3]
        treat.to_csv(path_or_buf=b.replace(
            ".bdg", "_subsIN.bdg"), sep="\t", index=False, header=None)


# 8. Standarize Samples


def standarizeBed(bed):

    df = pd.read_csv(bed, sep="\t", header=None)
    df.iloc[:, 3:] = df.iloc[:, 3:].apply(
        lambda x: (x-x.mean())/x.std(), axis=0)
    df.to_csv(path_or_buf=bed.replace(
        ".bdg", "_STZ.bdg"), sep="\t", index=False, header=None)


rpm_beds = [b for b in os.listdir(wd) if b.endswith("_subsIN.bdg")]

for bed in rpm_beds:
    standarizeBed(bed)

# 9. Create IGV tracks for the bdg files.

file_list = [b for b in os.listdir() if b.endswith(".bdg")]


def calltoTDF(infile):

    outfile = infile+".tdf"
    cmd = ("~/Programs/IGV_2.4.10/IGVTools/igvtools toTDF "
           f"{infile} {outfile} "
           "~/Programs/IGV_2.4.10/Custom_Genomes/PlasmoDB-41_Pfalciparum3D7.genome")

    sp.call(cmd, shell=True)


for file in file_list:
    calltoTDF(file)
