import sys
import SNPANN_upgrade_3

if __name__ == "__main__":
    #read input from stdin
    line = sys.stdin.readline()
    while line:
        line = line.strip()
        info_tuple = SNPANN_upgrade_3.SNP_ANN_of_Gene_Structure(line)
        ## output  tuple to stdout with tab as delimiter
        sys.stdout.write("\t".join(info_tuple) + "\n")

        line = sys.stdin.readline()