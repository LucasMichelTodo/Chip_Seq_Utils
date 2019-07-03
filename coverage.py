import os
import subprocess as sp
import pandas as pd
import numpy as np

# Functions


def removeDuplicates(bam):
    """ Call GATK remove duplicates function on a bam file. """

    i = indir+bam
    o = wd+bam.replace(".bam", "_noDup.bam")
    m = wd+bam.replace(".bam", "_metrics.txt")

    cmd = ("java -jar $PICARD MarkDuplicates "
           "REMOVE_DUPLICATES=true I={} O={} M={}") .format(i, o, m)

    sp.call(cmd, shell=True)


def indexBam(bam):
    """ Index a bamfile using samtools. """

    cmd = "samtools index {}" .format(bam)
    sp.call(cmd, shell=True)


def coverage(bam, out):
    """ Calculate per base coverage of a paired-end BAM file. """

    cmd = "bedtools genomecov -ibam {} -pc -d > {}" .format(bam, out)
    sp.call(cmd, shell=True)


def properBED(prebed):
    """ Reformat the output of bedtools genomecov to BED. """

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


def normalizeRPM(bed):
    """ Divide coverage in each position by number of reads. 
        Requires nreads, a dictinoray with bams as entries and 
        number of reads as values. """

    with open(bed, "r+") as infile:
        with open(bed.replace(".bdg", "_RPM.bdg"), "w+") as outfile:
            for line in infile:
                linelist = line.strip().split("\t")
                raw_val = int(linelist[3])
                genlen = int(linelist[2]) - int(linelist[1])

                bamname = bed.replace("_depth.bdg", ".bam")
                rpm = raw_val/(nreads[bamname]/2000000)

                outfile.write("\t".join(linelist[0:3])+"\t"+str(rpm)+"\n")


def divideByInput(chip_bed, input_bed, pseudo=0.1, log=True):
    """ Divide chip_bed signal by chip_input signal.
        Add a pseudocount to both chip and input.
        Optionally log the result (log2). """

    chip = pd.read_csv(chip_bed, sep="\t", header=None)
    ctl = pd.read_csv(input_bed, sep="\t", header=None)

    chip.iloc[:, 3] = (chip.iloc[:, 3]+pseudo)/(ctl.iloc[:, 3]+pseudo)

    if log:
        chip.iloc[:, 3] = np.log2(chip.iloc[:, 3])

    chip.to_csv(path_or_buf=chip_bed.replace(
        ".bdg", "_normIN.bdg"), sep="\t", index=False, header=None)


def substractInput(chip_bed, input_bed,):

    treat = pd.read_csv(chip_bed, sep="\t", header=None)
    ctl = pd.read_csv(input_bed, sep="\t", header=None)

    treat.iloc[:, 3] = treat.iloc[:, 3]-ctl.iloc[:, 3]
    treat.to_csv(path_or_buf=chip_bed.replace(
        ".bdg", "_subsIN.bdg"), sep="\t", index=False, header=None)


def standarizeBed(bed):

    df = pd.read_csv(bed, sep="\t", header=None)
    df.iloc[:, 3:] = df.iloc[:, 3:].apply(
        lambda x: (x-x.mean())/x.std(), axis=0)
    df.to_csv(path_or_buf=bed.replace(
        ".bdg", "_STZ.bdg"), sep="\t", index=False, header=None)


def calltoTDF(infile):
    """ Input a .bdg file (bed) and create an IGV readable .tdf using igvtools. """

    outfile = infile+".tdf"
    cmd = ("~/Programs/IGV_2.4.10/IGVTools/igvtools toTDF "
           f"{infile} {outfile} "
           "~/Programs/IGV_2.4.10/Custom_Genomes/PlasmoDB-41_Pfalciparum3D7.genome")

    sp.call(cmd, shell=True)


# Workflow
# All files must be in the same directory.
# All files must follow the naming convention: name_in/me/ac_...


indir = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/Renamed/RPM/"
wd = "/media/lucas/Disc4T/Projects/Cristina_ChIP_All/Coverage/Renamed/RPM/"
os.chdir(wd)

bams = [b for b in os.listdir(indir) if b.endswith("q5.bam")]

for bam in bams:
    removeDuplicates(bam)

bams = [b for b in os.listdir(wd) if b.endswith("_noDup.bam")]

for bam in bams:
    indexBam(bam)

nreads = {}
for b in bams:

    cmd = "samtools view -c {}" .format(b)
    nreads[b] = int(sp.check_output(cmd, shell=True))

for bam in bams:
    coverage(bam, bam.replace(".bam", "_depth.prebed"))

prebeds = [b for b in os.listdir(wd) if b.endswith("_depth.prebed")]

for pre in prebeds:
    properBED(pre)

beds = [b for b in os.listdir(wd) if b.endswith("_depth.bdg")]

for bed in beds:
    normalizeRPM(bed)

rpm_beds = [b for b in os.listdir(wd) if b.endswith("_RPM.bdg")]

for b in rpm_beds:
    mark = b.split("_")[1]
    prfx = b.split("_")[0]
    if mark == "me" or mark == "ac":
        treat_file = b
        ctl_file = [x for x in rpm_beds
                    if x.split("_")[0] == prfx and
                    x.split("_")[1] == "in"][0]

rpm_beds = [b for b in os.listdir(wd) if b.endswith("_subsIN.bdg")]

for bed in rpm_beds:
    standarizeBed(bed)

file_list = [b for b in os.listdir() if b.endswith(".bdg")]

for file in file_list:
    calltoTDF(file)
