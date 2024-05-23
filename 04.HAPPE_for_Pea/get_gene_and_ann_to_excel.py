import sys
import argparse
import os
import logging
gene_ann = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/gene.func.expression.pl"
python3 = "/public/home/fengcong/anaconda2/envs/py3/bin/python"
to_excel_pre = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/csv_to_Three-line_table_prepare.py"
to_excel = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/csv_to_Three-line_table.py"
gff_file = "/vol1/agis/chengshifeng_group/fengcong/30sample_popTest/04.snp_anno/IWGSC_v1.1_HC_20170706.gff3"

logging.basicConfig(level = logging.DEBUG,format = '[%(asctime)s] - %(name)s - %(levelname)s - %(message)s',stream = sys.stderr)
logger = logging.getLogger("Gene and Annotation")

def check_ret(ret,step):
    if  ret:
        logger.error("%s -- retrun code:%d\n"%(step,ret))
        exit(-1)

if __name__ == "__main__":
    cmdparser = argparse.ArgumentParser(description="get gene and annotation to excel/fengcong@caas.cn" )
    cmdparser.add_argument("-r","--region", dest="region",type=str, required=True,
                            help="region format , eg. chr1A:1-100.")
    cmdparser.add_argument("-g","--gff", dest="gff",type=str, required=False,default=gff_file , 
                            help="gff file" )
    cmdparser.add_argument("-o","--outputprefix", dest="outputprefix",type=str, required=True , 
                            help="output prefix." )

    args = cmdparser.parse_args()
    region = args.region
    gff = args.gff
    oprefix = args.outputprefix

    chromosome = region.split(":")[0]
    chr_start = int(region.split(":")[1].split("-")[0])
    chr_end = int(region.split(":")[1].split("-")[1])


    output_dir = os.path.dirname(oprefix)
    # print(output_dir)
    #create output dir
    if not os.path.exists(output_dir) and output_dir:
        os.makedirs(output_dir)

    ret = os.system("""
    awk '{if($3=="gene" && $1 == "%s" && $4 >= %d && $4 <= %d)print}'  %s | sed 's:.*\\tID=\\([^;]*\\);.*:\\1:g' > %s.gene.list
    
    """%(chromosome,chr_start,chr_end,gff, oprefix )
    )

    check_ret(ret,"get gene list.")


    ret = os.system("""
    perl %s %s.gene.list  | sed 's:,:;:g' | sed 's:\t:,:g' > %s.csv
    """%(gene_ann, oprefix , \
        oprefix ) 
    )

    check_ret(ret,"get gene ann.")

    ret = os.system("""
    %s %s %s.csv 1 1 "Table, Gene annotation and expression in %s."  > %s.csv.pre
    """%(python3,to_excel_pre ,oprefix, \
        region,oprefix) 
    )

    check_ret(ret,"to excel pre.")

    ret = os.system("""
    %s %s -o %s.geneann.xlsx %s.csv.pre 
    """%(python3, to_excel ,oprefix, oprefix) 
    )

    check_ret(ret,"to excel.")

    ##remove tmpfile
    os.remove("%s.gene.list"%(oprefix))
    os.remove("%s.csv"%(oprefix))
    os.remove("%s.csv.pre"%(oprefix))




