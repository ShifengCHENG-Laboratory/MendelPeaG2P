#!/usr/bin/env python3
# encoding: utf-8
"""
@version: v3.0
@author: Huangzejian
@project: common
@file: blocks_connection.py
@time: 2021/7/12 10:20
usage: python3 blocks_connection.py chr7B.ldblock.blocks.det  chr7B.SNP.vcf.gz 1000 25 0.98 /vol3/agis/chengshifeng/huangzejian/output/
"""

import random
from subprocess import call
import sys
import os
import gzip
import numpy as np
import json
import argparse
from multiprocessing.pool import Pool
from functools import partial

bcftools_cmd = "bcftools view -r {0} {1}"
plink_cmd = "plink --vcf {0} -r2 dprime with-freqs --allow-extra-chr --ld-window-r2 0 --out {1} --ld-window 1000000 --ld-window-kb 200000"
count = 0


def is_not_none(s):
    if s != "|":
        return s and s.strip()


def get_block_file_data(path):
    """
    从det文件中获取数据信息
    """
    dt = {}
    with open(path, 'r') as f:
        block_lines = f.readlines()[1:]
        for i in range(len(block_lines)):
            block_line_pre = block_lines[i].strip().split(" ")
            block_line = list(filter(is_not_none, block_line_pre))
            block_line = block_lines[i].strip().split("\t") if len(block_line) == 1 else block_line
            block_title = block_line[0] + "_" + block_line[1] + "_" + block_line[2]
            block_start = int(block_line[1])
            block_end = int(block_line[2])
            block_size = int(float(block_line[3]) * 1000)
            block_num_snps = int(block_line[4])
            block_snps = block_line[-1]
            dt[block_title] = [i, block_start, block_end, block_size, block_num_snps, block_snps]

    # sort: 根据i的值
    res = sorted(dt.items(), key=lambda k: k[1])
    return res


def get_window_value(snps_v, ix_lst):
    return [snps_v[i] for i in ix_lst]


def get_snps_id(window, start, end, snps_v, block_size):
    # 全部snp位点都选择
    return snps_v


def complete_snp_points(chr_name, points):
    return [chr_name + "_" + str(i) for i in points]


def handle_adjacent_block(item):
    # 对一条block的信息进行处理
    title_item = item[0]
    size_item = int(item[1][3])
    snps_item = item[1][-1]
    snps_item_value = [int(i.strip().split("_")[1]) for i in snps_item.split("|")]
    start = int(item[1][1])
    end = int(item[1][2])
    return [size_item, snps_item, snps_item_value, start, end, title_item]


def blocks_calc_points_select(binsize: int, data: list):
    points_data = {}
    chr_name = data[0][0].split("_")[0]
    for i in range(len(data)):
        if i < len(data) - 1:
            # [size_item, snps_item, snps_item_value, start, end, title_item]
            chr_a = handle_adjacent_block(item=data[i])
            chr_b = handle_adjacent_block(item=data[i + 1])

            # get_snps_id(window, start, end, snps_v, block_size)
            snps_a_id = get_snps_id(window=binsize, start=chr_a[3], end=chr_a[4], snps_v=chr_a[2], block_size=chr_a[0])
            snps_b_id = get_snps_id(window=binsize, start=chr_b[3], end=chr_b[4], snps_v=chr_b[2], block_size=chr_b[0])

            adjacent_block_select_points = []
            adjacent_block_select_points.extend(snps_a_id)
            adjacent_block_select_points.extend(snps_b_id)

            adjacent_block_select_points.sort()

            adjacent_block_key = str(chr_a[-1]) + "-vs-" + str(chr_b[-1])
            points_data[adjacent_block_key] = [
                complete_snp_points(chr_name, adjacent_block_select_points),
                [adjacent_block_select_points[0], adjacent_block_select_points[-1]], snps_a_id, snps_b_id,
                len(chr_a[2]),
                len(chr_b[2])
            ]
        else:
            pass
    return points_data


def execute_bcftools_cmd(_bcf_out_file, _arg1, _arg2, _block1, _block2):
    # 执行bcftools命令，并且过滤bcf生成的不在选取的中的位点
    with gzip.open(_bcf_out_file, 'wt') as gf:
        with os.popen(bcftools_cmd.format(_arg1, _arg2), 'r') as p:
            # bcftools_lines = p.readlines()
            for bcf_line in p:
                if bcf_line.startswith(chromosome):
                    bcf_line_lst = bcf_line.strip().split("\t")
                    snp_id_v = int(bcf_line_lst[2].split("_")[1])
                    if snp_id_v in _block1 or snp_id_v in _block2:
                        gf.write(bcf_line)
                    else:
                        pass
                else:
                    gf.write(bcf_line)


def statistics_plink_info(data_lst):
    quantile10 = np.percentile(data_lst, 10, interpolation='lower')
    lower = np.percentile(data_lst, 25, interpolation='lower')
    median = np.median(data_lst)
    upper = np.percentile(data_lst, 75, interpolation='lower')
    quantile90 = np.percentile(data_lst, 90, interpolation='lower')
    mean = np.mean(data_lst)
    return [quantile10, lower, median, upper, quantile90, mean]

def multi_handle_plink_line(handle_line):
    handle_line = handle_line.strip()
    # if not handle_line.strip().startswith(chromosome):
    if not handle_line.startswith("CHR_A"):
        line_lst_pre = handle_line.split(" ")
        line_lst = list(filter(is_not_none, line_lst_pre))
        bp_a = int(line_lst[1])
        bp_b = int(line_lst[5])
        if bp_a in block1 and bp_b in block1:
            pass
        elif bp_a in block2 and bp_b in block2:
            pass
        else:
            dp = float(line_lst[-1])
            r2 = float(line_lst[-2])
            return [dp, r2]
    else:
        return None


def process_file_line(_plink_out_file, line):
    results = []
    with open(_plink_out_file, 'r') as p1:
        for line in p1:
            result = multi_handle_plink_line(line)
            results.append(result)
    return results

def execute_plink_cmd(_plink_filter_out_file, _plink_out_file, _block1, _block2, _block1_snps_len, _block2_snps_len, _proportion, _plink_out_prefix):
    # 执行plink命令，并且过滤链接的位点同时存在在同一条block的情况。统计剩余位点的DP值。
    _dp_lst = []
    _r2_lst = []
    block_ls = _plink_out_prefix.split("-vs-")
    block1_len = int(block_ls[0].split("_")[2]) - int(block_ls[0].split("_")[1])
    block2_len = int(block_ls[1].split("_")[2]) - int(block_ls[1].split("_")[1])
    pos = int((int(block_ls[1].split("_")[1]) + int(block_ls[0].split("_")[2])) / 2)

    # M1:
    # with open(_plink_out_file, 'r') as p1:
    #     plink_lines = p1.readlines()
    #     with Pool(processes=30) as pool:
    #         rets = pool.map(multi_handle_plink_line, plink_lines)

    # M2:
    # 使用functools.partial将_plink_out_file作为参数预设传递给process_file_line
    partial_process_file_line = partial(process_file_line, _plink_out_file)

    with Pool(processes=30) as pool:
        rets = pool.map(partial_process_file_line, range(30))

    for ii in rets[0]:
        if ii is not None and ii:
            _dp_lst.append(ii[0])
            _r2_lst.append(ii[1])

    # [quantile10, lower, median, upper, quantile90, mean]
    if _dp_lst and _r2_lst:
        dp_ret = statistics_plink_info(_dp_lst)
        r2_ret = statistics_plink_info(_r2_lst)
        dp_lst_quantile = np.percentile(_dp_lst, _proportion, interpolation='lower')
        if dp_lst_quantile > confidence:
            connect = "y"
        else:
            connect = "n"

        ret = [
            pos, block1_len, _block1_snps_len,
            len(_block1), block2_len, _block2_snps_len,
            len(_block2), r2_ret[0], r2_ret[1], r2_ret[2], r2_ret[3], r2_ret[4], r2_ret[5], dp_ret[0], dp_ret[1],
            dp_ret[2], dp_ret[3], dp_ret[4], dp_ret[5], connect
        ]
        return ret
    else:
        return None


def block_confidence_judgment(prop, confide, dp_values):
    gt = [True for i in dp_values if float(i) > confide]
    if float(len(gt) / len(dp_values)) >= prop:
        return True
    else:
        return False


def get_2items_block_info(block_1, block_2):
    # 处理相邻的两个block，在进行两个block进行相连的时候
    conn_block_title = block_1[0] + "-vs-" + block_2[0]
    orderi = block_1[1][0]
    bp1 = block_1[1][1]
    bp2 = block_2[1][2]
    kb = block_2[1][2] - block_1[1][1] + 1
    nsnps = block_1[1][4] + block_2[1][4]
    snps = block_1[1][-1] + "|" + block_2[1][-1]
    return [conn_block_title, [orderi, bp1, bp2, kb, nsnps, snps]]


def get_item_block_info(block_1):
    # 处理一个block，在进行两个block进行相连的时候
    conn_block_title = block_1[0]
    orderi = block_1[1][0]
    bp1 = block_1[1][1]
    bp2 = block_1[1][2]
    kb = block_1[1][3]
    nsnps = block_1[1][4]
    snps = block_1[1][-1]
    return [conn_block_title, [orderi, bp1, bp2, kb, nsnps, snps]]


def blocks_merge(connect_block, origin_data):
    # 根据连续标志位判断两个block是否应该连续的合并下去；循环所有的block，将相邻的两个block连接，
    # 判断是否满足连接条件，若满足，进行连接，并且设置连续连接标志位为真，且将当前连接数据存放在连续标志位列表中；
    # 如果满足条件，判断连续标志位列表是否有数据，有，则取出数据进行数据合并，之后将连接数据存放在连续标志位列表中；无，则存放当前数据
    # 若不满足，不进行连接，并且设置连续连接标志位为假；
    connected_block_dt = []
    # origin_data = [(block_title, [i, block_start, block_end, block_size, block_num_snps, block_snps])]
    continuous_flag = False
    continuous_flag_dt = []
    for i in range(len(origin_data)):
        if i < len(origin_data) - 1:
            block_1 = origin_data[i]
            block_2 = origin_data[i + 1]
            prediction_blocks = block_1[0] + "-vs-" + block_2[0]
            if prediction_blocks in connect_block:
                two_block_info = get_2items_block_info(block_1, block_2)
                if continuous_flag:
                    if continuous_flag_dt:
                        block_conn_pre = continuous_flag_dt[0]
                        pre_blocks = block_conn_pre[0].split("-vs-")[-1]
                        current_blocks = block_1[0]
                        if pre_blocks == current_blocks:
                            _two_block_info = get_2items_block_info(block_conn_pre, block_2)
                            # _nsnps = block_2[1][4] + block_conn_pre[1][4]
                            continuous_flag_dt = [(_two_block_info[0], [
                                _two_block_info[1][0], _two_block_info[1][1], _two_block_info[1][2],
                                _two_block_info[1][3], _two_block_info[1][4], _two_block_info[1][5]
                            ])]
                            continuous_flag = True
                        else:
                            continuous_flag = False
                            connected_block_dt.append(continuous_flag_dt[0])
                            continuous_flag_dt = []
                else:
                    if continuous_flag_dt:
                        continuous_flag_dt.append((two_block_info[0], [
                            two_block_info[1][0], two_block_info[1][1], two_block_info[1][2], two_block_info[1][3],
                            two_block_info[1][4], two_block_info[1][5]
                        ]))
                        continuous_flag = True
                    else:
                        continuous_flag_dt = [(two_block_info[0], [
                            two_block_info[1][0], two_block_info[1][1], two_block_info[1][2], two_block_info[1][3],
                            two_block_info[1][4], two_block_info[1][5]
                        ])]
                        continuous_flag = True

            else:
                one_block_info = get_item_block_info(block_1)
                if continuous_flag_dt:
                    _block_conn_pre = continuous_flag_dt[0]
                    connected_block_dt.append(_block_conn_pre)
                    continuous_flag_dt = []
                else:
                    connected_block_dt.append((one_block_info[0], [
                        one_block_info[1][0], one_block_info[1][1], one_block_info[1][2], one_block_info[1][3],
                        one_block_info[1][4], one_block_info[1][5]
                    ]))
                continuous_flag = False
        else:
            block_1 = origin_data[-1]
            one_block_info = get_item_block_info(block_1)
            if continuous_flag_dt:
                _block_conn_pre = continuous_flag_dt[0]
                connected_block_dt.append(_block_conn_pre)
                _pre_blocks = _block_conn_pre[0].split("-vs-")[-1]
                if _pre_blocks != one_block_info[0]:
                    connected_block_dt.append((one_block_info[0], [
                        one_block_info[1][0], one_block_info[1][1], one_block_info[1][2], one_block_info[1][3],
                        one_block_info[1][4], one_block_info[1][5]
                    ]))
            else:
                connected_block_dt.append((one_block_info[0], [
                    one_block_info[1][0], one_block_info[1][1], one_block_info[1][2], one_block_info[1][3],
                    one_block_info[1][4], one_block_info[1][5]
                ]))

    return connected_block_dt


def blocks_conn(origin_data, chr_out_mean_file, final_out_det_file, _confidence):
    # ld_dt = {}
    connect_block = []
    with open(chr_out_mean_file, 'r') as ldf:
        ldf_liens = ldf.readlines()
        for ldf_line in ldf_liens[1:]:
            line_lst = ldf_line.strip().split("\t")
            # line_lst = list(filter(is_not_none, line_lst_pre))
            group = line_lst[0]
            _connect = line_lst[-1]
            # 判断两个block是否满足90%的点DP值大于98%
            if _connect == 'y':
                connect_block.append(group)

    # 进行相邻block的连接
    _chromosome = origin_data[0][0].split("_")[0]
    ret_data = blocks_merge(connect_block, origin_data)
    res = sorted(ret_data, key=lambda k: k[1])
    with open(final_out_det_file, 'w') as ff:
        ff.write("CHR" + "\t" + "BP1" + "\t" + "BP2" + "\t" + "KB" + "\t" + "NSNPS" + "\t" + "SNPS" + "\n")
        for item in res:
            ff.write(
                str(_chromosome) + "\t" + str(item[1][1]) + "\t" + str(item[1][2]) + "\t" +
                str(float(item[1][3]) / 1000) + "\t" + str(item[1][4]) + "\t" + str(item[1][5]) + "\n")


if __name__ == "__main__":
    # args
    import time
    start = time.time()
    parser = argparse.ArgumentParser(description="Implement the connection between the blocks.")
    parser.add_argument("-B", "--block_file", type=str, required=True, help="The blocks file for the chromosome.")
    parser.add_argument("-v", "--vcf_file", type=str, required=True, help="The vcf file for the chromosome(.gz).")
    parser.add_argument(
        "-b",
        "--bin_size",
        type=int,
        required=True,
        help="Calculate the number of points required to select for each block.")
    parser.add_argument(
        "-q",
        "--quantile",
        type=int,
        required=True,
        help="Select the confidence quantile for the two block sites of the plink connection.")
    parser.add_argument(
        "-c",
        "--confidence",
        type=float,
        required=True,
        help="The confidence selected by the two block sites for the plink connection.")
    parser.add_argument("-o", "--out_dir", type=str, required=True, help="Output folder.")

    args = parser.parse_args()

    block_file = args.block_file
    vcf_gz_file = args.vcf_file
    bin_size = args.bin_size
    quantile = args.quantile
    confidence = args.confidence
    output = args.out_dir

    if not os.path.exists(output):
        os.system("mkdir -p {0}".format(output))

    block_data = get_block_file_data(path=block_file)
    chr_points_data = blocks_calc_points_select(binsize=bin_size, data=block_data)

    chromosome = list(chr_points_data.keys())[0].split("-vs-")[0].split("_")[0]
    total_out_mean_file = os.path.join(output, chromosome + ".total.mean")
    # total_out_file = os.path.join(output, chromosome + ".total.ld")
    final_block_conn_file = os.path.join(output, chromosome + ".ldblock.blocks.final.det")

    os.system(
        "echo '{0}	{1}	{2}	{3}	{4}	{5}	{6}	{7}	{8}	{9}	{10}	{11}	{12}	{13}	{14}	{15}	{16}	{17}	{18}	{19}	{20}' > {21}"
        .format("group", "pos", "block1_len", "block1_snps_num", "block1_snps_select", "block2_len", "block2_snps_num",
                "block2_snps_select", "quantile10_R2", "quantile25_R2", "quantile50_R2", "quantile75_R2",
                "quantile90_R2", "mean_R2", "quantile10_DP", "quantile25_DP", "quantile50_DP", "quantile75_DP",
                "quantile90_DP", "mean_DP", "connect", total_out_mean_file))
    # os.system("> {0}".format(total_out_file))

    for key, value in chr_points_data.items():
        # chromosome = key.split("-vs-")[0].split("_")[0]
        start_item = value[1][0]
        end_item = value[1][1]
        block1 = value[2]
        block2 = value[3]
        block1_snps_len = value[4]
        block2_snps_len = value[5]
        arg1 = chromosome + ":" + str(start_item) + "-" + str(end_item)
        arg2 = vcf_gz_file

        bcf_out_file = os.path.join(output, chromosome + ".vcf.gz")
        execute_bcftools_cmd(_bcf_out_file=bcf_out_file, _arg1=arg1, _arg2=arg2, _block1=block1, _block2=block2)

        plink_filter_out_file = os.path.join(output, key + ".filter.ld")
        plink_out_file = os.path.join(output, key + ".ld")
        plink_out_prefix = os.path.join(output, key)

        call(plink_cmd.format(bcf_out_file, plink_out_prefix), shell=True)
        # pf = open(plink_filter_out_file, 'wt')
        plink_lst = execute_plink_cmd(
            _plink_filter_out_file=plink_filter_out_file,
            _plink_out_file=plink_out_file,
            _block1=block1,
            _block2=block2,
            _block1_snps_len=block1_snps_len,
            _block2_snps_len=block2_snps_len,
            _proportion=quantile,
            _plink_out_prefix=key)

        if plink_lst:
            # os.system("echo '{0} {1} {2}' >> {3}".format(key, arr_mean, arr_var, total_out_mean_file))
            # [pos, block1_len, _block1_snps_len, len(_block1), block2_len, _block1_snps_len, len(_block2), dp_ret[0], dp_ret[1], dp_ret[2],
            #   dp_ret[3], dp_ret[4], dp_ret[5], r2_ret[0], r2_ret[1], r2_ret[2], r2_ret[3], r2_ret[4], r2_ret[5], connect]
            os.system(
                "echo '{0}	{1}	{2}	{3}	{4}	{5}	{6}	{7}	{8}	{9}	{10}	{11}	{12}	{13}	{14}	{15}	{16}	{17}	{18}	{19}	{20}' >> {21}"
                .format(key, plink_lst[0], plink_lst[1], plink_lst[2], plink_lst[3], plink_lst[4], plink_lst[5],
                        plink_lst[6], plink_lst[7], plink_lst[8], plink_lst[9], plink_lst[10], plink_lst[11],
                        plink_lst[12], plink_lst[13], plink_lst[14], plink_lst[15], plink_lst[16], plink_lst[17],
                        plink_lst[18], plink_lst[19], total_out_mean_file))

        # os.system("cat {0} >> {1}".format(plink_out_file, total_out_file))
        os.system("rm -f {0} {1} {2}".format(plink_out_file, plink_out_prefix + "*", bcf_out_file))
        # pf.close()
    sys.stdout.write("******chromosome ld blocks filter done!******\n")

    sys.stdout.write("******start blocks connection******\n")
    blocks_conn(
        origin_data=block_data,
        chr_out_mean_file=total_out_mean_file,
        final_out_det_file=final_block_conn_file,
        _confidence=confidence)
    sys.stdout.write("******blocks connection done!******\n")
    sys.stdout.write("******all program done!******\n")
    os.system("touch {0}".format(chromosome + ".done"))
    end = time.time()
    sys.stdout.write("******Total time consuming: {0}******\n".format(end - start))
