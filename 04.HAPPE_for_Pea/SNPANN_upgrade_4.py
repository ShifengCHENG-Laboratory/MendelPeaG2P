# coding:utf-8
# 注释选择规则：首先按照feature_id末尾的数字进行分组（依次分为转录本1，转录本2...等n个组）；在各组内依据注释优先级排名，注释优先级见anno_sort_rule)，两次排序后取第一个为该snp唯一注释。
import re
import pandas as pd


def SNP_ANN_of_Gene_Structure(INFO_field):
    """
    @INFO_field:vcf INFO field, which contains SNP information. should be containd ANN field.
    @return: a annotation of SNP in gene structure.
    @returntype: str
    return value one of as follows:
    None object : error
    tuple object : featureid，annotation，c.xx,p.xx
    ...
    """
    
    # 判断是否含有ANN字段
    if 'ANN=' in INFO_field:
        pass
    else:
        return None
    
    return_dict = {
        "3_prime_UTR_variant":
            "3_prime_UTR_variant",
        "5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant":
            "start_gain",
        "5_prime_UTR_variant":
            "5_prime_UTR_variant",
        "5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|splice_region_variant":
            "start_gain",
        "3_prime_UTR_variant|splice_region_variant":
            "3_prime_UTR_variant",
        "splice_region_variant|3_prime_UTR_variant":
            "3_prime_UTR_variant",
        "splice_region_variant|5_prime_UTR_variant":
            "5_prime_UTR_variant",
        "5_prime_UTR_variant|splice_region_variant":
            "5_prime_UTR_variant",
        "downstream_gene_variant":
            "downstream_gene_variant",
        "intergenic_region":
            "intergenic_region",
        "intron_variant":
            "intron_variant",
        "splice_acceptor_variant&intron_variant":
            "splice_acceptor_variant",
        "splice_donor_variant&intron_variant":
            "splice_donor_variant",
        "splice_region_variant&intron_variant":
            "splice_region_variant",
        "splice_acceptor_variant&splice_donor_variant&intron_variant":
            "splice_acceptor_variant&splice_donor_variant",
        "splice_acceptor_variant&splice_region_variant&intron_variant":
            "splice_acceptor_variant",
        "splice_donor_variant&splice_region_variant&intron_variant":
            "splice_donor_variant",
        "splice_region_variant&stop_retained_variant":
            "stop_retained_variant",
        "start_lost&splice_region_variant":
            "start_lost",
        "stop_lost&splice_region_variant":
            "stop_lost",
        "splice_region_variant&synonymous_variant":
            "splice_region_variant",
        "missense_variant":
            "missense_variant",
        "missense_variant&splice_region_variant":
            "missense_variant",
        "start_lost":
            "start_lost",
        "stop_gained":
            "stop_gained",
        "stop_gained&splice_region_variant":
            "stop_gained",
        "stop_lost":
            "stop_lost",
        "stop_retained_variant":
            "stop_retained_variant",
        "synonymous_variant":
            "synonymous_variant",
        "upstream_gene_variant":
            "upstream_gene_variant"
    }
    df_mapping = pd.DataFrame({
        'anno_sort_rule': [
            'stop_gained', 'stop_gained&splice_region_variant', 'stop_lost',
            'stop_lost&splice_region_variant', 'start_lost',
            'start_lost&splice_region_variant',
            'splice_acceptor_variant&intron_variant',
            'splice_donor_variant&intron_variant', 'missense_variant',
            'missense_variant&splice_region_variant',
            'splice_acceptor_variant&splice_region_variant&intron_variant',
            'splice_donor_variant&splice_region_variant&intron_variant',
            'splice_region_variant&synonymous_variant',
            'stop_retained_variant',
            'splice_region_variant&stop_retained_variant',
            'synonymous_variant', '5_prime_UTR_variant|splice_region_variant',
            '5_prime_UTR_variant', '3_prime_UTR_variant',
            '3_prime_UTR_variant|splice_region_variant',
            '5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant',
            '5_prime_UTR_premature_start_codon_gain_variant|5_prime_UTR_variant|splice_region_variant',
            'splice_acceptor_variant&splice_donor_variant&intron_variant',
            'splice_region_variant&intron_variant', 'intron_variant',
            'upstream_gene_variant', 'downstream_gene_variant',
            'intergenic_region'
        ]
    })
    num_df_mapping = pd.DataFrame({
        'num_sort_rule':
            ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0', 'D']
    })
    sort_mapping = df_mapping.reset_index().set_index('anno_sort_rule')
    num_sort_mapping = num_df_mapping.reset_index().set_index('num_sort_rule')
    inf_line = INFO_field.split(';')  # vcf文件第八列的注释信息以；分隔
    anno_dic = {}  # {id1:[[ann1,ann2],(c.a>t),(p.his>tyr)],[],...}
    for anno_block in inf_line:
        if re.match('ANN', anno_block):  # 从头匹配ANN，详细注释从这里开始
            ANN_block = anno_block[4:].split(',')  # 多个注释以，分隔
            #    feature_id_list = []  # [id1,id2,...]
            for anno_inf in ANN_block:
                anno_inf = anno_inf.split("|")  # 第一个|后即alle的类别信息,第7列为feature id
                annotation = anno_inf[1]
                impact = anno_inf[2]
                feature_type = anno_inf[5]
                feature_id = anno_inf[8]
                cxx = anno_inf[11]
                pxx = anno_inf[12]
                if feature_id not in anno_dic:
                    anno_dic[feature_id] = {'Annotation': [], 'c.xx': set(), 'p.xx': set(), 'impact': set()}
                if annotation == 'intergenic_region':
                    anno_dic[feature_id]['Annotation'].append('intergenic_region')
                    anno_dic[feature_id]['c.xx'].add(cxx)
                    anno_dic[feature_id]['p.xx'].add(pxx)
                    anno_dic[feature_id]['impact'].add(impact)
                elif annotation == "intragenic_variant":
                    anno_dic[feature_id]['Annotation'].append("intragenic_variant")
                    anno_dic[feature_id]['c.xx'].add(cxx)
                    anno_dic[feature_id]['p.xx'].add(pxx)
                    anno_dic[feature_id]['impact'].add(impact)
                elif feature_type == "transcript":
                    anno_dic[feature_id]['Annotation'].append(annotation)
                    anno_dic[feature_id]['c.xx'].add(cxx)
                    anno_dic[feature_id]['p.xx'].add(pxx)
                    anno_dic[feature_id]['impact'].add(impact)
                else:
                    anno_dic.pop(feature_id, None)
    # 把每个snp的注释存入dataframe，最后依据dataframe的不同列进行排序。
    final_all_anno = [
        [_id, '|'.join(anno_dic[_id]['Annotation']), '|'.join(anno_dic[_id]['c.xx']), '|'.join(anno_dic[_id]['p.xx']),
         '|'.join(anno_dic[_id]['impact'])] for _id in anno_dic
    ]
    
    snp_ann_df = pd.DataFrame(final_all_anno, columns=['Feature_id', 'Annotation', 'c.xx', 'p.xx', 'impact'])
    snp_ann_df['second_sort'] = snp_ann_df['Annotation'].map(sort_mapping['index'])
    snp_ann_df['last_num'] = snp_ann_df['Feature_id'].str[-1]
    snp_ann_df['first_sort'] = snp_ann_df['last_num'].map(num_sort_mapping['index'])
    
    sort_ann_df = snp_ann_df.sort_values(by=['first_sort', 'second_sort'], ascending=[True, True]).reset_index(drop=True)

    # 以元组形式输出（featureid，annotation，c.xx,p.xx)
    if not sort_ann_df.empty:
        out_anno = (
            sort_ann_df.loc[0, 'Feature_id'],
            return_dict.get(sort_ann_df.loc[0, 'Annotation'], sort_ann_df.loc[0, 'Annotation']),
            sort_ann_df.loc[0, 'c.xx'],
            sort_ann_df.loc[0, 'p.xx'],
            sort_ann_df.loc[0, 'impact']
        )
    else:
        out_anno = (None, None, None, None, None)
    return out_anno
