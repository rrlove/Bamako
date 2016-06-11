##this is a test of Git with this IDE
#!/usr/bin/env python

Usage = """
vcf-consensus.py - Becca Love, 2 June 2016

Take a VCF file and a reference genome, and construct a consensus sequence from a particular sample with sensitivity to read counts. For example, a SNP called as a heterozygote with reference/alternate read counts of 6/4 is more likely to be incorporated into the reference sequence than a SNP called as 10/3.

Importantly, this is a stochastic process and the same consensus sequence will not be constructed each time from the same sample and reference genome.

This was designed for use with VCFs from pooled data, where read counts may be much higher.

Optionally, the SNPs may be filtered during this process by passing the argument -q with a minimum site quality, corresponding to field 6 of the VCF.

Currently, this script only works on biallelic sites.

Usage:
        vcf-consensus.py -v vcf -o outprefix -s sample -r reference -q min_SNP_quality
        
Depends:
        Biopython https://biopython.org/DIST/docs/tutorial/Tutorial.html
        vcf/pysam https://pysam.readthedocs.io/en/latest/api.html

"""

##to-do: implement support for multiple samples in a population

import argparse
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord
import random
import unittest
import vcf

parser = argparse.ArgumentParser(description="Generate a consensus sequence from one or more samples in a VCF file, taking into account minor allele frequencies.")

parser.add_argument('-v','--vcf', action="store", dest="vcfname", type=str, help='VCF file used to generate consensus sequence. Required.', required=True)
parser.add_argument('-o','--output', action="store", dest="outbase", type=str, help='Prefix for output files. If no prefix is provided, the sample name specified is used.')##note below, default is sample name
parser.add_argument('-s','--sample', action="store", dest="sample", type=str, help='1 sample name for which to generate a consensus. Required.', required=True)
parser.add_argument('-r','--reference', action="store", dest="reference", type=str, help='The reference sequence against which the VCF file was generated. Required.', required=True)
parser.add_argument('-q','--quality',action="store",dest="quality", type=float, help='Minimum quality for considering a locus.')
parser.add_argument('-x', '--verbose', action="store", dest="verbosity", help='Run the program with progress messages.')

args = parser.parse_args()

if not args.outbase:
    args.outbase = args.sample

def analyze_locus(record,sample):
    
    altFreq = 0
    cumAlts = 0
    cumReads = 0
    
    ##generate a random number to which to compare the alternate allele frequency
    testBound = random.uniform(0.0,1.0)##should I use scipy.random.rand() here again?
    
    if not record.genotype(sample)["AD"]:
        print "Warning: missing genotype in " + sample + " at " + record.CHROM + ":" + str(record.POS) + "; using reference allele"
        return record.REF
    
    elif record.genotype(sample)["AD"]:
        
        if len(record.genotype(sample)["AD"]) == 2:
            ##for a two-allele locus, the # of reads for the alternate allele will always be 2nd in the list generated by ["AD"].
            cumAlts += float(record.genotype(sample)["AD"][1])    ##get # reads supporting alternate allele
            cumReads += float(sum(record.genotype(sample)["AD"]))    ##get total # reads
    
        elif len(record.genotype(sample)["AD"] > 2):    ##only biallelic loci are supported for now
            ##to-do: accommodate multiallelic loci
            print "Skipping multi-allelic locus at " + record.CHROM + ":" + str(record.POS)
        
        altFreq = float(cumAlts) / float(cumReads)
    
        if altFreq > testBound:
            return record.ALT[0]
        elif altFreq <= testBound:
            return record.REF

vcf_reader = vcf.Reader(open(args.vcfname,'r'))

reference = SeqIO.index(args.reference, "fasta")

try:
    ##make sure that the reference genome and the supplied VCF file have the same contigs
    assert set(vcf_reader.contigs.keys()) == set(reference.keys())
except:
    print "Error: contigs of reference and VCF do not match"
    raise SystemExit

##set a minimum quality to 30 x # samples if one is not supplied
if not args.quality:
    args.quality = 30 * len(vcf_reader.samples)

##make the Seq objects for the modified contigs that will get written, and a dictionary to hold them
outContigs = {}

for contig in vcf_reader.contigs.keys():
    outContigs[contig] = SeqRecord.SeqRecord("",id=contig,name=contig)

##keep track of where in the reference we last stopped
stoppedAt = 0

for record in vcf_reader:
    
    if record.POS < stoppedAt:##to-do: make sure this catches transition between contigs
        stoppedAt = 0
    
    ##add everything up to the position before the current locus
    outContigs[record.CHROM] += reference[record.CHROM].seq[stoppedAt:(record.POS-1)]
    ##make sure the reference alleles match in the reference genome and the VCF
    assert record.REF == reference[record.CHROM].seq[record.POS-1].upper(), "Error: reference and VCF do not match"    ##check this for off-by-one error
    
    ##ignore filtered sites and sites not passing the quality threshold
    if record.QUAL >= args.quality and not record.FILTER:
    
        outContigs[record.CHROM].seq += analyze_locus(record,args.sample)
        stoppedAt = record.POS+len(record.REF)-1##check this for off-by-one error,account for indels!
    
    if args.verbosity:
        print "Finished processing " + record.CHROM + ":" + str(record.POS)
        
outFile = args.outbase + ".fasta"

for contigs in outContigs.keys():
    if outContigs[contig].seq:
        SeqIO.write(outContigs[contig], outFile, "fasta")
    print contig