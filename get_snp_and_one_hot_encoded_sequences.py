import argparse
import numpy as np
from ucscgenome import Genome

def oneHotEncodeSequence(sequence):
    oneHotDimension = (len(sequence), 4)
    dnaAlphabet = {"A":0, "G":1, "C":2, "T":3}    
    one_hot_encoded_sequence = np.zeros(oneHotDimension, dtype=np.int)
    for i, nucleotide in enumerate(sequence):
        if nucleotide.upper() in dnaAlphabet:
            index = dnaAlphabet[nucleotide.upper()]
            one_hot_encoded_sequence[i][index] = 1
    return one_hot_encoded_sequence


def getSequence(allele,chrom,position,left,right,genome):
    alleleLength = len(allele)
    deductable = "right"
    for i in range(alleleLength-1):
        if deductable=="right":
            right-=1
            deductable="left"
        elif deductable=="left":
            left-=1
            deductable="right"
    
    left_sequence = genome[chrom][position-left:position]
    right_sequence = genome[chrom][position+alleleLength:position+alleleLength+right]
    
    sequence = left_sequence.lower()+allele.lower()+right_sequence.lower()
    return sequence


def getOneHotEncodedSequences(bedInput, xOutput, xAlternateOutput, leftWindow, rightWindow, genomeObject, pValueCutoff, snpDataFile):
    if snpDataFile:
        outf = open(snpDataFile, 'w')
    finalReferenceSequences = []
    finalAlternateSequences = []
    with open(bedInput, 'r') as f:
        for line in f:
            curInterval = line.strip().split()
            chrom = curInterval[0]
            start = int(curInterval[1])
            end = int(curInterval[2])
            rsid = curInterval[3]
            a1 = curInterval[4]
            a2 = curInterval[5]
            gwas_p = float(curInterval[6])
            gwas_z = float(curInterval[7])
            if gwas_p < pValueCutoff:
                if snpDataFile:
                    outf.write("\t".join(curInterval))
                    outf.write("\n")
                referenceSequence  = getSequence(a1,chrom,start,leftWindow,rightWindow,genomeObject)
                alternateSequence = getSequence(a2,chrom,start,leftWindow,rightWindow,genomeObject)
                assert(len(referenceSequence)==leftWindow+rightWindow+1)
                assert(len(alternateSequence)==leftWindow+rightWindow+1)
                encodedReferenceSequence = oneHotEncodeSequence(referenceSequence)
                encodedAlternateSequence = oneHotEncodeSequence(alternateSequence)
                finalReferenceSequences.append(encodedReferenceSequence)
                finalAlternateSequences.append(encodedAlternateSequence)
            
    finalReferenceSequences = np.stack(finalReferenceSequences, axis=0)
    finalAlternateSequences = np.stack(finalAlternateSequences, axis=0)

    print finalReferenceSequences.shape
    print finalAlternateSequences.shape
    np.save(xOutput, finalReferenceSequences)
    np.save(xAlternateOutput, finalAlternateSequences)
    if snpDataFile:
        outf.close()
    
    
if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Construct numpy arrays containing one-hot encoded sequences for given snp bed file along with window size (BE CAREFUL IF WINDOWS ARE TOO BIG, SNP METADATA MAY NOT MATCH THE NUMPY ARRAYS)")
    parser.add_argument('-i', '--bed-like-input', help='bed file containing coordinates for snps, alleles, p and z scores', required=True)
    parser.add_argument('-xo', '--x-output', help = 'output numpy array of one hot encoded sequences', required=True)
    parser.add_argument('-xao', '--x-alternate-output', help = 'output numpy array of one hot encoded sequences for alternate alleles', required=True)
    parser.add_argument('-l', '--left-window', type=int, help = 'how many nucleotides to pad at the left', required=True)
    parser.add_argument('-r', '--right-window', type=int, help = 'how many nucleotides to pad at the right', required=True)    
    parser.add_argument('-g', '--genome-name', help='name of the genome', required=True)
    parser.add_argument('-d', '--genome-dir', help='local path to genomes (MUST BE in 2-bit format)', required=True)
    parser.add_argument('-p', '--p-value-cutoff', help='p value cutoff', type=float, required=False, default=1.0)
    parser.add_argument('-s', '--snp-data-file', help='output snp bed file with allele information', required=False)
    
    args = parser.parse_args()
    bedInput = args.bed_like_input
    xOutput = args.x_output
    xAlternateOutput = args.x_alternate_output
    leftWindow = args.left_window
    rightWindow = args.right_window
    genomeName = args.genome_name
    genomeDir = args.genome_dir
    pValueCutoff = args.p_value_cutoff
    snpDataFile = args.snp_data_file
    genomeObject = Genome(genomeName, cache_dir=genomeDir, use_web=False)
    getOneHotEncodedSequences(bedInput, xOutput, xAlternateOutput, leftWindow, rightWindow, genomeObject, pValueCutoff, snpDataFile)