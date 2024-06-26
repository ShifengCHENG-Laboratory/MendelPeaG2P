# coding:utf-8
# 注释选择规则：首先按照feature_id末尾的数字进行分组（依次分为转录本1，转录本2...等n个组）；在各组内依据注释优先级排名，注释优先级见anno_sort_rule)，两次排序后取第一个为该Indel唯一注释。indel注释仅分为exon和intron，downstream以及upstream
import re
import pandas as pd


def INDEL_ANN_of_Gene_Structure(INFO_field):
    '''
    @INFO_field:vcf INFO field, which contains SNP information. should be containd ANN field.
    @return: a annotation of SNP in gene structure.
    @returntype: str
    return value one of as follows:
    None object : error
    ...
    '''

    # 判断是否含有ANN字段
    if 'ANN=' in INFO_field:
        pass
    else:
        return None

    return_dict = {
        "bidirectional_gene_fusion":
        "bidirectional_gene_fusion",
        "downstream_gene_variant":
        "downstream_gene_variant",
        "downstream_gene_variant|frameshift_variant&stop_gained":
        "downstream_gene_variant",
        "3_prime_UTR_truncation&exon_loss_variant&splice_region_variant":
        "exon_variant",
        "3_prime_UTR_variant":
        "3_prime_UTR_variant",
        "3_prime_UTR_variant|frameshift_variant&stop_lost":
        "exon_variant",
        "3_prime_UTR_variant|frameshift_variant&stop_lost&splice_region_variant":
        "exon_variant",
        "3_prime_UTR_variant|non_coding_transcript_variant|splice_region_variant":
        "exon_variant",
        "3_prime_UTR_variant|splice_region_variant":
        "splice_region_variant",
        "3_prime_UTR_variant|stop_lost&conservative_inframe_deletion":
        "exon_variant",
        "3_prime_UTR_variant|stop_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost":
        "exon_variant",
        "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost&splice_region_variant":
        "exon_variant",
        "5_prime_UTR_variant":
        "5_prime_UTR_variant",
        "5_prime_UTR_variant|frameshift_variant&start_lost":
        "exon_variant",
        "5_prime_UTR_variant|frameshift_variant&start_lost&splice_region_variant":
        "exon_variant",
        "5_prime_UTR_variant|splice_region_variant":
        "splice_region_variant",
        "5_prime_UTR_variant|start_lost&conservative_inframe_deletion":
        "exon_variant",
        "5_prime_UTR_variant|start_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "conservative_inframe_deletion":
        "exon_variant",
        "conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "conservative_inframe_insertion":
        "exon_variant",
        "conservative_inframe_insertion&splice_region_variant":
        "exon_variant",
        "disruptive_inframe_deletion":
        "exon_variant",
        "disruptive_inframe_deletion&splice_region_variant":
        "exon_variant",
        "disruptive_inframe_insertion":
        "exon_variant",
        "disruptive_inframe_insertion&splice_region_variant":
        "exon_variant",
        "exon_loss_variant|frameshift_variant":
        "exon_variant",
        "exon_loss_variant&splice_region_variant":
        "exon_variant",
        "frameshift_variant":
        "exon_variant",
        "frameshift_variant&start_lost":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant":
        "exon_variant",
        "frameshift_variant&stop_gained":
        "exon_variant",
        "frameshift_variant&stop_lost":
        "exon_variant",
        "frameshift_variant&stop_lost&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant":
        "exon_variant",
        "frameshift_variant&stop_lost&splice_region_variant":
        "exon_variant",
        "start_lost&conservative_inframe_deletion":
        "exon_variant",
        "start_lost&conservative_inframe_insertion":
        "exon_variant",
        "start_lost&conservative_inframe_insertion&splice_region_variant":
        "exon_variant",
        "start_lost&disruptive_inframe_deletion":
        "exon_variant",
        "start_lost&disruptive_inframe_insertion":
        "exon_variant",
        "start_lost&disruptive_inframe_insertion&splice_region_variant":
        "exon_variant",
        "stop_gained&conservative_inframe_insertion":
        "exon_variant",
        "stop_gained&conservative_inframe_insertion&splice_region_variant":
        "exon_variant",
        "stop_gained&disruptive_inframe_deletion":
        "exon_variant",
        "stop_gained&disruptive_inframe_deletion&splice_region_variant":
        "exon_variant",
        "stop_gained&disruptive_inframe_insertion":
        "exon_variant",
        "stop_gained&disruptive_inframe_insertion&splice_region_variant":
        "exon_variant",
        "stop_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "stop_lost&conservative_inframe_insertion":
        "exon_variant",
        "stop_lost&conservative_inframe_insertion&splice_region_variant":
        "exon_variant",
        "stop_lost&disruptive_inframe_deletion":
        "exon_variant",
        "stop_lost&disruptive_inframe_deletion&splice_region_variant":
        "exon_variant",
        "stop_lost&disruptive_inframe_insertion":
        "exon_variant",
        "stop_lost&disruptive_inframe_insertion&splice_region_variant":
        "exon_variant",
        "stop_lost&conservative_inframe_deletion":
        "exon_variant",
        "frameshift_variant&stop_gained&splice_region_variant":
        "exon_variant",
        "frameshift_variant&splice_region_variant":
        "exon_variant",
        "3_prime_UTR_variant|downstream_gene_variant|non_coding_transcript_variant|splice_region_variant":
        "exon_variant",
        "3_prime_UTR_truncation&exon_loss_variant&splice_region_variant|downstream_gene_variant|non_coding_transcript_variant":
        "exon_variant",
        "3_prime_UTR_variant|non_coding_transcript_variant|splice_region_variant&downstream_gene_variant":
        "exon_variant",
        "conservative_inframe_deletion&splice_region_variant|downstream_gene_variant":
        "exon_variant",
        "downstream_gene_variant|exon_loss_variant&splice_acceptor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "downstream_gene_variant|frameshift_variant&splice_region_variant":
        "exon_variant",
        "downstream_gene_variant|frameshift_variant&stop_lost":
        "exon_variant",
        "downstream_gene_variant|frameshift_variant&stop_lost&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant":
        "exon_variant",
        "downstream_gene_variant|frameshift_variant&stop_lost&splice_region_variant":
        "exon_variant",
        "downstream_gene_variant|stop_lost&3_prime_UTR_truncation&exon_loss_variant&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "downstream_gene_variant|stop_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "exon_loss_variant&splice_acceptor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "exon_loss_variant&splice_donor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&stop_lost&splice_acceptor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "stop_lost&splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_donor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant|splice_donor_variant&5_prime_UTR_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant|splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&stop_lost&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&stop_lost&splice_region_variant|splice_acceptor_variant&3_prime_UTR_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&stop_lost&splice_region_variant|splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant":
        "exon_variant",
        "splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|stop_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "splice_donor_variant&5_prime_UTR_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant":
        "exon_variant",
        "stop_gained&splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "stop_gained&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "stop_lost&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&conservative_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_acceptor_variant&splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_acceptor_variant&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "splice_donor_variant&disruptive_inframe_deletion&intron_variant":
        "exon_variant",
        "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant":
        "exon_variant",
        "downstream_gene_variant|non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&downstream_gene_variant&intron_variant":
        "exon_variant",
        "exon_loss_variant&splice_donor_variant&splice_region_variant&intron_variant|upstream_gene_variant":
        "exon_variant",
        "exon_loss_variant&splice_region_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|upstream_gene_variant":
        "exon_variant",
        "frameshift_variant|upstream_gene_variant":
        "exon_variant",
        "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost|upstream_gene_variant":
        "exon_variant",
        "5_prime_UTR_truncation&exon_loss_variant|non_coding_transcript_variant|upstream_gene_variant":
        "exon_variant",
        "5_prime_UTR_truncation&exon_loss_variant|start_lost&conservative_inframe_deletion|upstream_gene_variant":
        "exon_variant",
        "5_prime_UTR_variant|non_coding_transcript_variant|upstream_gene_variant":
        "5_prime_UTR_variant",
        "5_prime_UTR_variant|upstream_gene_variant":
        "5_prime_UTR_variant",
        "conservative_inframe_deletion|upstream_gene_variant":
        "exon_variant",
        "conservative_inframe_insertion|upstream_gene_variant":
        "exon_variant",
        "exon_loss_variant|frameshift_variant&start_lost|upstream_gene_variant":
        "exon_variant",
        "exon_loss_variant&splice_region_variant|upstream_gene_variant":
        "exon_variant",
        "start_lost&conservative_inframe_deletion|upstream_gene_variant":
        "exon_variant",
        "frameshift_variant&start_lost&splice_region_variant|upstream_gene_variant":
        "exon_variant",
        "frameshift_variant&start_lost|upstream_gene_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&5_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_acceptor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&3_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant|upstream_gene_variant":
        "exon_variant",
        "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant&upstream_gene_variant":
        "exon_variant",
        "gene_fusion":
        "gene_fusion",
        "intergenic_region":
        "intergenic_region",
        "intragenic_variant":
        "intragenic_variant",
        "intron_variant":
        "intron_variant",
        "splice_acceptor_variant&intron_variant":
        "splice_region_variant",
        "splice_acceptor_variant&splice_donor_variant&intron_variant":
        "splice_region_variant",
        "splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant":
        "splice_region_variant",
        "splice_acceptor_variant&splice_region_variant&intron_variant":
        "splice_region_variant",
        "splice_donor_variant&intron_variant":
        "splice_region_variant",
        "splice_donor_variant&splice_region_variant&intron_variant":
        "splice_region_variant",
        "splice_region_variant&intron_variant":
        "splice_region_variant",
        "regulatory_region_variant":
        "regulatory_region_variant",
        "upstream_gene_variant":
        "upstream_variant"
    }
    df_mapping = pd.DataFrame({
        'anno_sort_rule': [
            "bidirectional_gene_fusion", "gene_fusion",
            "frameshift_variant&stop_lost&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "downstream_gene_variant|frameshift_variant&stop_lost&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant",
            "frameshift_variant&stop_lost&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant",
            "frameshift_variant&start_lost&splice_region_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "frameshift_variant&start_lost&splice_region_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&intron_variant",
            "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost&splice_region_variant",
            "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost|upstream_gene_variant",
            "5_prime_UTR_truncation&exon_loss_variant|frameshift_variant&start_lost",
            "exon_loss_variant|frameshift_variant&start_lost|upstream_gene_variant",
            "exon_loss_variant|frameshift_variant",
            "stop_lost&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
            "downstream_gene_variant|stop_lost&3_prime_UTR_truncation&exon_loss_variant&conservative_inframe_deletion&splice_region_variant",
            "splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant",
            "5_prime_UTR_truncation&exon_loss_variant|start_lost&conservative_inframe_deletion|upstream_gene_variant",
            "exon_loss_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "downstream_gene_variant|exon_loss_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
            "exon_loss_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant",
            "downstream_gene_variant|non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&downstream_gene_variant&intron_variant",
            "exon_loss_variant&splice_region_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|upstream_gene_variant",
            "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant|upstream_gene_variant",
            "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant|splice_region_variant&upstream_gene_variant",
            "non_coding_transcript_variant|splice_donor_variant&3_prime_UTR_truncation&exon_loss_variant&splice_region_variant&intron_variant",
            "exon_loss_variant&splice_donor_variant&splice_region_variant&intron_variant|upstream_gene_variant",
            "exon_loss_variant&splice_donor_variant&splice_region_variant&intron_variant",
            "5_prime_UTR_truncation&exon_loss_variant|non_coding_transcript_variant|upstream_gene_variant",
            "3_prime_UTR_truncation&exon_loss_variant&splice_region_variant|downstream_gene_variant|non_coding_transcript_variant",
            "exon_loss_variant&splice_region_variant|upstream_gene_variant",
            "exon_loss_variant&splice_region_variant",
            "3_prime_UTR_truncation&exon_loss_variant&splice_region_variant",
            "frameshift_variant&stop_gained&splice_region_variant",
            "frameshift_variant&stop_gained",
            "frameshift_variant&stop_lost&splice_region_variant|splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant",
            "frameshift_variant&stop_lost&splice_region_variant|splice_acceptor_variant&3_prime_UTR_variant&intron_variant",
            "frameshift_variant&stop_lost&splice_acceptor_variant&splice_region_variant&intron_variant",
            "3_prime_UTR_variant|frameshift_variant&stop_lost&splice_region_variant",
            "3_prime_UTR_variant|frameshift_variant&stop_lost",
            "downstream_gene_variant|frameshift_variant&stop_lost&splice_region_variant",
            "downstream_gene_variant|frameshift_variant&stop_lost",
            "frameshift_variant&stop_lost&splice_region_variant",
            "frameshift_variant&stop_lost",
            "frameshift_variant&start_lost&splice_region_variant|splice_donor_variant&5_prime_UTR_variant&intron_variant",
            "frameshift_variant&start_lost&splice_region_variant|splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant",
            "frameshift_variant&start_lost&splice_donor_variant&splice_region_variant&intron_variant",
            "5_prime_UTR_variant|frameshift_variant&start_lost&splice_region_variant",
            "frameshift_variant&start_lost&splice_region_variant|upstream_gene_variant",
            "frameshift_variant&start_lost&splice_region_variant",
            "5_prime_UTR_variant|frameshift_variant&start_lost",
            "frameshift_variant&start_lost|upstream_gene_variant",
            "frameshift_variant&start_lost",
            "frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant",
            "frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant",
            "frameshift_variant&splice_donor_variant&splice_region_variant&intron_variant",
            "downstream_gene_variant|frameshift_variant&splice_region_variant",
            "frameshift_variant&splice_region_variant",
            "frameshift_variant|upstream_gene_variant", "frameshift_variant",
            "stop_gained&splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "stop_gained&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "stop_gained&disruptive_inframe_insertion",
            "stop_gained&disruptive_inframe_insertion&splice_region_variant",
            "stop_gained&conservative_inframe_insertion&splice_region_variant",
            "stop_gained&conservative_inframe_insertion",
            "stop_gained&disruptive_inframe_deletion",
            "stop_gained&disruptive_inframe_deletion&splice_region_variant",
            "stop_lost&splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|stop_lost&conservative_inframe_deletion&splice_region_variant",
            "stop_lost&disruptive_inframe_insertion&splice_region_variant",
            "stop_lost&disruptive_inframe_insertion",
            "stop_lost&conservative_inframe_insertion&splice_region_variant",
            "stop_lost&conservative_inframe_insertion",
            "stop_lost&disruptive_inframe_deletion&splice_region_variant",
            "stop_lost&disruptive_inframe_deletion",
            "3_prime_UTR_variant|stop_lost&conservative_inframe_deletion&splice_region_variant",
            "stop_lost&conservative_inframe_deletion&splice_region_variant",
            "3_prime_UTR_variant|stop_lost&conservative_inframe_deletion",
            "downstream_gene_variant|stop_lost&conservative_inframe_deletion&splice_region_variant",
            "stop_lost&conservative_inframe_deletion",
            "splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant",
            "splice_donor_variant&5_prime_UTR_variant&intron_variant|start_lost&conservative_inframe_deletion&splice_region_variant",
            "start_lost&disruptive_inframe_insertion&splice_region_variant",
            "start_lost&disruptive_inframe_insertion",
            "start_lost&conservative_inframe_insertion&splice_region_variant",
            "start_lost&conservative_inframe_insertion",
            "start_lost&disruptive_inframe_deletion",
            "5_prime_UTR_variant|start_lost&conservative_inframe_deletion&splice_region_variant",
            "5_prime_UTR_variant|start_lost&conservative_inframe_deletion",
            "start_lost&conservative_inframe_deletion|upstream_gene_variant",
            "start_lost&conservative_inframe_deletion",
            "splice_acceptor_variant&splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "splice_acceptor_variant&splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant",
            "splice_acceptor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
            "splice_acceptor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&5_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_acceptor_variant&3_prime_UTR_variant&intron_variant|splice_region_variant",
            "splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&intron_variant",
            "splice_donor_variant&disruptive_inframe_deletion&intron_variant",
            "splice_donor_variant&conservative_inframe_deletion&splice_region_variant&intron_variant",
            "non_coding_transcript_variant|splice_donor_variant&splice_region_variant&5_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_donor_variant&splice_region_variant&3_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_donor_variant&5_prime_UTR_variant&intron_variant|splice_region_variant",
            "non_coding_transcript_variant|splice_donor_variant&3_prime_UTR_variant&intron_variant|splice_region_variant",
            "disruptive_inframe_insertion&splice_region_variant",
            "disruptive_inframe_insertion",
            "conservative_inframe_insertion&splice_region_variant",
            "conservative_inframe_insertion|upstream_gene_variant",
            "conservative_inframe_insertion",
            "disruptive_inframe_deletion&splice_region_variant",
            "disruptive_inframe_deletion",
            "conservative_inframe_deletion&splice_region_variant|downstream_gene_variant",
            "conservative_inframe_deletion&splice_region_variant",
            "conservative_inframe_deletion|upstream_gene_variant",
            "conservative_inframe_deletion",
            "5_prime_UTR_variant|splice_region_variant",
            "3_prime_UTR_variant|non_coding_transcript_variant|splice_region_variant",
            "3_prime_UTR_variant|splice_region_variant",
            "3_prime_UTR_variant|downstream_gene_variant|non_coding_transcript_variant|splice_region_variant",
            "3_prime_UTR_variant|non_coding_transcript_variant|splice_region_variant&downstream_gene_variant",
            "5_prime_UTR_variant|non_coding_transcript_variant|upstream_gene_variant",
            "5_prime_UTR_variant|upstream_gene_variant", "5_prime_UTR_variant",
            "3_prime_UTR_variant", "intragenic_variant",
            "splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant",
            "splice_acceptor_variant&splice_donor_variant&intron_variant",
            "splice_acceptor_variant&splice_region_variant&intron_variant",
            "splice_acceptor_variant&intron_variant",
            "splice_donor_variant&splice_region_variant&intron_variant",
            "splice_donor_variant&intron_variant",
            "splice_region_variant&intron_variant", "intron_variant",
            "upstream_gene_variant",
            "downstream_gene_variant|frameshift_variant&stop_gained",
            "downstream_gene_variant", "intergenic_region",
            "regulatory_region_variant"
        ]
    })
    num_df_mapping = pd.DataFrame({
        'num_sort_rule':
        ['1', '0', 'D', '2', '3', '4', '5', '6', '7', '8', '9']
    })
    sort_mapping = df_mapping.reset_index().set_index('anno_sort_rule')
    num_sort_mapping = num_df_mapping.reset_index().set_index('num_sort_rule')
    inf_line = INFO_field.split(';')  # vcf文件第八列的注释信息以；分隔
    for anno_block in inf_line:
        if re.match('ANN', anno_block):  # 从头匹配ANN，详细注释从这里开始
            ANN_block = anno_block[4:].split(',')  # 多个注释以，分隔
            anno_dic = {}  # {id1:[[ann1,ann2],(c.a>t),(p.his>tyr)],[],...}
            #    feature_id_list = []  # [id1,id2,...]
            for anno_inf in ANN_block:
                anno_inf = anno_inf.split(
                    "|")  # 第一个|后即alle的类别信息,第7列为feature id
                annotation = anno_inf[1]
                feature_id = anno_inf[8]
                cxx = anno_inf[10]
                pxx = anno_inf[11]
                feature_type = anno_inf[5]
                if feature_id not in anno_dic:
                    anno_dic[feature_id] = []
                    anno_dic[feature_id].append([])
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id].append(set())
                if annotation == 'intergenic_region':
                    anno_dic.pop(feature_id)
                    feature_id=anno_inf[6]
                    anno_dic[feature_id] = []
                    anno_dic[feature_id].append([])
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id][0].append('intergenic_region')
                    anno_dic[feature_id][1].add(cxx)
                    anno_dic[feature_id][2].add(pxx)
                elif annotation == 'bidirectional_gene_fusion':
                    # 这个变异的id那列不全，应该重新赋值feature_id才能知道受影响的两个基因分别是哪个。
                    anno_dic.pop(feature_id)
                    feature_id = anno_inf[4]
                    anno_dic[feature_id] = []
                    anno_dic[feature_id].append([])
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id][0].append('bidirectional_gene_fusion')
                    anno_dic[feature_id][1].add(cxx)
                    anno_dic[feature_id][2].add(pxx)
                elif annotation == 'gene_fusion':
                    anno_dic.pop(feature_id)
                    feature_id = anno_inf[4]
                    anno_dic[feature_id] = []
                    anno_dic[feature_id].append([])
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id].append(set())
                    anno_dic[feature_id][0].append('gene_fusion')
                    anno_dic[feature_id][1].add(cxx)
                    anno_dic[feature_id][2].add(pxx)
                elif annotation == "intragenic_variant":
                    anno_dic[feature_id][0].append('intragenic_variant')
                    anno_dic[feature_id][1].add(cxx)
                    anno_dic[feature_id][2].add(pxx)
                # elif feature_id.startswith("Traes") and feature_id[-2] == '.':
                elif feature_type == "transcript":
                    anno_dic[feature_id][0].append(annotation)
                    anno_dic[feature_id][1].add(cxx)
                    anno_dic[feature_id][2].add(pxx)
                else:
                    #这一步必须保留，否则后面排序会被regulatory及对应id干扰。
                    anno_dic.pop(feature_id)
    # 把每个snp的注释存入dataframe，最后依据dataframe的不同列进行排序。
    snp_ann_df = pd.DataFrame()
    final_all_anno = []
    for id in anno_dic:  # 将每个feature_id对应的所有注释转化为1个
#        print(id)
        anno_dic[id][0].sort()
        anno_list = anno_dic[id][0]
        all_anno = '|'.join(anno_list)
        cxx_list = anno_dic[id][1]
        all_cxx = '|'.join(cxx_list)
        pxx_list = anno_dic[id][2]
        all_pxx = '|'.join(pxx_list)
        final_all_anno.append([id, all_anno, all_cxx,
                               all_pxx])  # 利用列表增加元素较dataframe直接增加快很多
        # anno_dic[id][0] = all_anno
        # anno_dic[id][1] = all_cxx
        # anno_dic[id][2] = all_pxx
        #anno_row = {
        #    'Feature_id': id,
        #    'Annotation': all_anno,
        #    'c.xx': all_cxx,
        #    'p.xx': all_pxx
    # }
    # snp_ann_df = snp_ann_df.append(anno_row, ignore_index=True) #循环加入非常耗时
    snp_ann_df = pd.DataFrame(
        final_all_anno, columns=['Feature_id', 'Annotation', 'c.xx', 'p.xx'])
    snp_ann_df['second_sort'] = snp_ann_df['Annotation'].map(
        sort_mapping['index'])
    snp_ann_df['contains_gene'] = snp_ann_df['Feature_id'].str.contains('gene')
    snp_ann_df['last_num'] = snp_ann_df['Feature_id'].map(lambda x: x[-1])
    snp_ann_df['first_sort'] = snp_ann_df['last_num'].map(
        num_sort_mapping['index'])

    # 先依据feature_id末尾的值排序，而按注释优先级排序
    sort_ann_df = snp_ann_df.sort_values(by=['contains_gene','first_sort', 'second_sort'],
                                         ascending=[True,True, True])
    sort_ann_df = sort_ann_df.reset_index(drop=True)
    out_anno = (sort_ann_df.loc[0, 'Feature_id'],
                return_dict[sort_ann_df.loc[0, 'Annotation']] if sort_ann_df.loc[0, 'Annotation'] in return_dict else sort_ann_df.loc[0, 'Annotation'],
                sort_ann_df.loc[0, 'c.xx'], sort_ann_df.loc[0, 'p.xx'])
    return out_anno  # 以元组形式输出（featureid，annotation，c.xx,p.xx)
