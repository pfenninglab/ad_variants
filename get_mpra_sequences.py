import pandas as pd
import argparse

import py2bit


def getUpdatedPaddings(allele,left,right):
    alleleLength = len(allele)
    
    deductable = "right"
    for i in range(alleleLength-1):
        if deductable=="right":
            right-=1
            deductable="left"
        elif deductable=="left":
            left-=1
            deductable="right"
    
    return left,right


def getSequencesAndPrint(dataFile, genome_object, mpraSeqLength, exptName, summitColumnPresent):
    if summitColumnNotPresent:
        data = pd.read_csv(dataFile,
                           sep='\t',
                           header=None,
                           names=["CHROM",
                                  "SNPSTART",
                                  "SNPEND",
                                  "RSID",
                                  "REF",
                                  "ALT",
                                  "PEAKCHROM",
                                  "PEAKSTART",
                                  "PEAKEND",
                                  "NAME",
                                  "SCORE",
                                  "_",
                                 ]
                          )
        data["SUMMIT"] = data["PEAKEND"]-data["PEAKSTART"]
        data["SUMMIT"] = data["SUMMIT"]//2
    else:
        data = pd.read_csv(dataFile,
                           sep='\t',
                           header=None,
                           names=["CHROM",
                                  "SNPSTART",
                                  "SNPEND",
                                  "RSID",
                                  "REF",
                                  "ALT",
                                  "PEAKCHROM",
                                  "PEAKSTART",
                                  "PEAKEND",
                                  "NAME",
                                  "SCORE",
                                  "STRAND",
                                  "SIGNAL",
                                  "P",
                                  "Q",
                                  "SUMMIT"
                                 ]
                          )        
    print("CHROM",
          "RSID",
          "SNPSTART",
          "SNPEND",
          "REF",
          "ALT",
          "EXPERIMENT",
          "PEAKSTART",
          "PEAKEND",
          "SUMMIT",
          "SUMMITZEROBASED",
          "SUMMITREFSEQ",
          "SUMMITALTSEQ",
          "SNPREFSEQ",
          "SNPALTSEQ",
          sep='\t'
     )
    for index,row in data.iterrows():
        chrom = row["CHROM"]
        ref = row["REF"]
        alt = row["ALT"]
        summitLocation = int(row["PEAKSTART"]+ row["SUMMIT"])
        mpraSeqStart = summitLocation - mpraSeqLength//2
        mpraSeqEnd = mpraSeqStart + mpraSeqLength
        snpStart = row["SNPSTART"]
        snpEnd = row["SNPEND"]
        if mpraSeqStart <= snpStart and snpEnd <= mpraSeqEnd and (snpStart+len(alt)) <= mpraSeqEnd:
            #making sure both the ref and alt alleles are within 227bp of summit
            (snpCenteredRefSequenceLeftWindow,
            snpCenteredRefSequenceRightWindow) = getUpdatedPaddings(ref,
                                                                    mpraSeqLength//2,
                                                                    mpraSeqLength//2
                                                                    )
            (snpCenteredAltSequenceLeftWindow,
            snpCenteredAltSequenceRightWindow) = getUpdatedPaddings(alt,
                                                                    mpraSeqLength//2,
                                                                    mpraSeqLength//2
                                                                    )
            snpCenteredRefSequenceLeft = genome_object.sequence(chrom,
                                                                snpStart-snpCenteredRefSequenceLeftWindow,
                                                                snpStart)
            snpCenteredRefSequenceRight = genome_object.sequence(chrom,
                                                                 snpStart+len(ref),
                                                                 snpStart+len(ref)+snpCenteredRefSequenceRightWindow
                                                                 )
            snpCenteredRefSequence = snpCenteredRefSequenceLeft.lower()+ref.lower()+snpCenteredRefSequenceRight.lower()



            snpCenteredAltSequenceLeft = genome_object.sequence(chrom,
                                                                snpStart-snpCenteredAltSequenceLeftWindow,
                                                                snpStart
                                                                )
            snpCenteredAltSequenceRight = genome_object.sequence(chrom,
                                                                 snpStart+len(ref),
                                                                 snpStart+len(ref)+snpCenteredAltSequenceRightWindow
                                                                 )
            snpCenteredAltSequence = snpCenteredAltSequenceLeft.lower()+alt.lower()+snpCenteredAltSequenceRight.lower()

            if snpStart < summitLocation:
                # case when snp is to the left of summit
                summitCenteredRefSequence = genome_object.sequence(chrom, snpEnd, mpraSeqEnd).lower()
                summitCenteredRefSequence = ref.lower() + summitCenteredRefSequence
                curRefSeqLen = len(summitCenteredRefSequence)
                if not curRefSeqLen == mpraSeqLength:
                    #edge case when the SNP is at the edge of the sequence
                    summitCenteredRefSequence = genome_object.sequence(chrom,
                                                                       snpStart+curRefSeqLen-mpraSeqLength,
                                                                       snpStart).lower() + summitCenteredRefSequence

                summitCenteredAltSequence = genome_object.sequence(chrom, snpEnd, mpraSeqEnd).lower()
                summitCenteredAltSequence = alt.lower()+summitCenteredAltSequence
                curAltSeqLen = len(summitCenteredAltSequence)
                if not curAltSeqLen == mpraSeqLength:
                    #edge case when the SNP is at the edge of the sequence
                    summitCenteredAltSequence = genome_object.sequence(chrom,
                                                                       snpStart+curAltSeqLen-mpraSeqLength,
                                                                       snpStart).lower() + summitCenteredAltSequence
            elif snpStart >= summitLocation:
                # case when snp is to the right of summit
                summitCenteredRefSequence = genome_object.sequence(chrom,mpraSeqStart,snpStart).lower()
                summitCenteredRefSequence = summitCenteredRefSequence + ref.lower()
                curRefSeqLen = len(summitCenteredRefSequence)
                if not curRefSeqLen == mpraSeqLength:
                    #edge case when the SNP is at the edge of the sequence
                    summitCenteredRefSequence = (summitCenteredRefSequence
                                                 +genome_object.sequence(chrom,
                                                                         snpEnd,
                                                                         snpEnd-curRefSeqLen+mpraSeqLength)).lower()

                summitCenteredAltSequence = genome_object.sequence(chrom,mpraSeqStart,snpStart).lower()
                summitCenteredAltSequence = summitCenteredAltSequence + alt.lower()
                curAltSeqLen = len(summitCenteredAltSequence)
                if not curAltSeqLen >= mpraSeqLength:
                    #edge case when the SNP is at the edge of the sequence
                    summitCenteredAltSequence = (summitCenteredAltSequence
                                                 +genome_object.sequence(chrom,
                                                                         snpEnd,
                                                                         snpEnd-curAltSeqLen+mpraSeqLength)).lower()



            assert(len(summitCenteredRefSequence)==mpraSeqLength)
            assert(len(summitCenteredAltSequence)==mpraSeqLength)
            assert(summitCenteredRefSequence.lower()==genome_object.sequence(chrom,
                                                                             mpraSeqStart,
                                                                             mpraSeqEnd).lower())

            print(row["CHROM"],
                  row["RSID"],
                  row["SNPSTART"],
                  row["SNPEND"],
                  row["REF"],
                  row["ALT"],
                  exptName,
                  row["PEAKSTART"],
                  row["PEAKEND"],
                  row["SUMMIT"],
                  summitLocation,
                  summitCenteredRefSequence,
                  summitCenteredAltSequence,
                  snpCenteredRefSequence,
                  snpCenteredAltSequence,
                  sep='\t'
                 )
                                             
                                             
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='get sequences for MPRA array', fromfile_prefix_chars='@')
    parser.add_argument('-n', '--name', help='name of cell type or experiment', required=True)
    parser.add_argument('-i', '--input', help='input text file (bedtools intersection of SNP file and peak file with -wa -wb)',
                        required=True)
    parser.add_argument('-s',
                        '--summit-column-not-present',
                        action='store_true',
                        help='include flag if summit column is not present in the intersection file. If true, middle of peak is treated as summit. False by default',
                        default=False)
    parser.add_argument('-l', '--seq-length', help='length of MPRA sequence', type=int, default=227)
    parser.add_argument('-g', '--genome', help='path to genome 2bit file', default="/home/eramamur/resources/genomes/hg19/hg19.2bit")                                          
    args = parser.parse_args()
    
    exptName = args.name
    dataFile = args.input                                             
    genome_object = py2bit.open(args.genome)
    mpraSeqLength = args.seq_length
    summitColumnNotPresent = args.summit_column_not_present                                             
    getSequencesAndPrint(dataFile, genome_object, mpraSeqLength, exptName, summitColumnNotPresent)
                                                 
    
