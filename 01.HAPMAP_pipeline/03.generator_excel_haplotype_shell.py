#!/usr/bin/env python3
# encoding: utf-8
"""
@version: v3.0
@author: Huangzejian
@project: common
@file: generator_excel_haplotype_sh.py
@time: 2021/9/26 10:20
usage: python3 blocks_connection/generator_excel_haplotype_sh.py
"""

import sys
import os
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format="[%(asctime)s] - [%(module)s - %(funcName)s - line:%(lineno)d] - %(levelname)s - %(message)s",
    stream=sys.stderr)
logger = logging.getLogger("generator_excel_haplotype_shell")

py3 = "/public/home/fengcong/anaconda2/envs/py3/bin/python"
# excel_py = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/excel_haplotype_update/excel_haplotype.py"
excel_py = "/public/agis/chengshifeng_group/fengcong/WGRS/software/Fc-code/excel_haplotype_update_only_clustering/excel_haplotype.py"
vcf_gz = "/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/03-ann/01-SNP/00-IC/{0}/{0}.snp.raw.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.InbreedingCoeff_retain.ann.vcf.gz"
inf_txt = "/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/002.SNP_recalling/01.SNP/08-hapmap/03.excel_Haplotype/zz_code/773_Inf_yr.txt"


def extract_info(path, out_dir):
    with open(path, 'r') as f:
        lines = f.readlines()

        # mkdir -p / 1000个block一个文件夹
        length = len(lines) - 1
        step = 1000
        dir_num = int(length / step) if length % step == 0 else int(length / step) + 1

        for i in range(1, dir_num + 1):
            start = (i - 1) * step + 1
            end = i * step
            if end > length:
                end = length
            if start == 1 and end >= 1000:
                dirs = "blocks_00001" + "_0" + str(end)
            elif start == 1 and end < 1000:
                dirs = "blocks_00001" + "_00" + str(end)
            elif start >= 1000 and end < 10000:
                dirs = "blocks_0" + str(start) + "_0" + str(end)
            else:
                dirs = "blocks_" + str(start) + "_" + str(end)
            block_path = os.path.join(out_dir, dirs)
            os.system("mkdir -p {0}".format(block_path))

            for i in range(start, end + 1):
                line = lines[i]
                line_lst = line.strip().split("\t")
                # logger.info(f"{line_lst}")
                chromosome = line_lst[0]
                region = line_lst[1] + "-" + line_lst[2]
                region_str = chromosome + ":" + region
                if i < 1000:
                    if len(str(i)) == 1:
                        block_name = "block0000" + str(i)
                    elif len(str(i)) == 2:
                        block_name = "block000" + str(i)
                    elif len(str(i)) == 3:
                        block_name = "block00" + str(i)
                elif i >= 1000 and i <= 10000:
                    block_name = "block0" + str(i)
                else:
                    block_name = "block" + str(i)
                out_name = chromosome + "_" + line_lst[1] + "_" + line_lst[2]
                out = block_name + "_" + out_name
                # functional
                # cmd_line = "{0} {1} -v {2} -i {3} -r {4} -F -o {5}".format(py3, excel_py, vcf_gz.format(chromosome), inf_txt, r_str, out)
                cmd_line = "{0} {1} -v {2} -i {3} -r {4} -o {5}".format(py3, excel_py, vcf_gz.format(chromosome), inf_txt, region_str, out)
                shell_path = os.path.join(block_path,  out_name + ".go.block.sh")
                os.system("echo '{0}' > {1}".format(cmd_line, shell_path))


if __name__ == "__main__":
    output_dir = "/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/09-hapmap/03.excel_haplotype/"
    in_base_path = "/vol3/agis/chengshifeng_group/shiyan/03-pea_WGS/003.SNP_recalling_773_V2/09-hapmap/02.LDblock_identify/output/{0}.ldblock.blocks.final.det"
    logger.info("start read infomations.")
    for i in range(5, 6):
        _chromosome = "chr" + str(i)
        in_file = in_base_path.format(_chromosome, _chromosome)
        chromosome_dir = os.path.join(output_dir, _chromosome)
        extract_info(in_file, chromosome_dir)
        logger.info("Chromosome {0} Shell script is generated.".format(_chromosome))
    logger.info("Done.")
