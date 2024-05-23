
    #!/usr/bin/sh
    source /public/home/fengcong/anaconda2/etc/profile.d/conda.sh
    conda init
    conda activate py3
    /public/home/fengcong/anaconda2/envs/py3/bin/python /public/agis/chengshifeng_group/fengcong/WGRS/software/gwas_pipline/merge_jpg_dirversion.py ./


    conda deactivate
    