import pandas as pd
import argparse
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
rcParams['svg.fonttype'] = 'none'
rcParams['font.size']=15




if __name__=="__main__":
    parser = argparse.ArgumentParser(description='get correlations between replicates', fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', help='input text file (bedtools subtract of klab unthresholded output and blacklist)', required=True)
    args = parser.parse_args()
    inFile = args.input

    rep1_rep2_data = pd.read_csv(inFile,
                                 sep="\t", compression="gzip",
                                 names=["CHR",
                                 "START",
                                 "END",
                                 "NAME",
                                 "SCORE",
                                 "STRAND",
                                 "SIGNAL",
                                 "P",
                                 "Q",
                                 "SUMMIT",
                                 "LOCALIDR",
                                 "GLOBALIDR",
                                 "REP1_START",
                                 "REP1_END",
                                 "REP1_SIGNAL",
                                 "REP1_SUMMIT",
                                 "REP2_START",
                                 "REP2_END",
                                 "REP2_SIGNAL",
                                 "REP2_SUMMIT"]
                                )
    print("\t".join(["Rep_pair",
                     "Pearson_all",
                     "Pearson_all_Pvalue",
                     "Spearman_all",
                     "Spearman_all_Pvalue",
                     "Pearson_chr4_validation",
                     "Pearson_chr4_validation_Pvalue",
                     "Spearman_chr4_validation",
                     "Spearman_chr4_validation_Pvalue",
                     "Pearson_train_chr",
                     "Pearson_train_chr_Pvalue",
                     "Spearman_train_chr",
                     "Spearman_train_chr_Pvalue"]))

    print(os.path.basename(inFile), sep="\t", end="\t")
    p_all = pearsonr(rep1_rep2_data["REP1_SIGNAL"], rep1_rep2_data["REP2_SIGNAL"])
    s_all = spearmanr(rep1_rep2_data["REP1_SIGNAL"], rep1_rep2_data["REP2_SIGNAL"])
    print(p_all[0], p_all[1], sep="\t", end="\t")
    print(s_all[0], s_all[1], sep="\t", end="\t")

    chr4_rep1_rep2_data = rep1_rep2_data.loc[rep1_rep2_data["CHR"]=="chr4"]
    p_4 = pearsonr(chr4_rep1_rep2_data["REP1_SIGNAL"], chr4_rep1_rep2_data["REP2_SIGNAL"])
    s_4 = spearmanr(chr4_rep1_rep2_data["REP1_SIGNAL"], chr4_rep1_rep2_data["REP2_SIGNAL"])
    print(p_4[0], p_4[1], sep="\t", end="\t")
    print(s_4[0], s_4[1], sep="\t", end="\t")


    chrs_training_rep1_rep2_data = rep1_rep2_data.loc[~rep1_rep2_data["CHR"].isin(["chr4", "chr8", "chr9"])]
    p_tr = pearsonr(chrs_training_rep1_rep2_data["REP1_SIGNAL"], chrs_training_rep1_rep2_data["REP2_SIGNAL"])
    s_tr = spearmanr(chrs_training_rep1_rep2_data["REP1_SIGNAL"], chrs_training_rep1_rep2_data["REP2_SIGNAL"])
    print(p_tr[0], p_tr[1], sep="\t", end="\t")
    print(s_tr[0], s_tr[1], sep="\t", end="\t")
    print()


    plt.scatter(chr4_rep1_rep2_data["REP1_SIGNAL"], chr4_rep1_rep2_data["REP2_SIGNAL"], s=5, alpha=0.5)
    plt.ylabel("DNase signal (Monocyte rep1)")
    plt.xlabel("DNase signal (Monocyte rep2)")
    plt.savefig("best_replicate_concordance.svg")
    plt.savefig("best_replicate_concordance.png")
    plt.close()
