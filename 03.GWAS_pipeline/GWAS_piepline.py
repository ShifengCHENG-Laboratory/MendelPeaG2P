#!/public/home/fengcong/anaconda2/bin/python
#--coding:utf-8--
#2020-5-14 15:39:29
#GWAS pipline
#fengcong@caas.cn
import sys
import gzip
import argparse
import os
from datetime import datetime
import subprocess


script_path = os.path.abspath(__file__)
manplot=os.path.dirname(script_path)+"/"+"manhattanplot_function.r"
qqplot=os.path.dirname(script_path)+"/"+"pretty_QQ_plot_func.r"
mergepy=os.path.dirname(script_path)+"/"+"merge_jpg_dirversion.py"
padj = os.path.dirname(script_path)+"/"+"padjust_and_r2.R"
pypath="/public/home/fengcong/anaconda2/envs/py3/bin/python"
rpath="/public/home/fengcong/anaconda2/envs/R/bin/Rscript"
emmaxpath="/public/agis/chengshifeng_group/fengcong/WGRS/software/emmax/"
farmcpujob=7

def go_farmcpupp(numfile,pheno,outdir):
    bedfile_abpath=numfile
    pheno_abpath=pheno

    if numfile.startswith("/"):
        #absolute path
        pass
    else:
        cwd=os.getcwd()
        bedfile_abpath = cwd+"/"+numfile
    
    if pheno.startswith("/"):
        #absolute path
        pass
    else:
        cwd=os.getcwd()
        pheno_abpath = cwd+"/"+pheno

    cmd="""
    ln -s %s*.txt %s;
    ln -s %s %s;
    """ %(bedfile_abpath,outdir,pheno_abpath,outdir)
    os.system(cmd)

    # generate script and run
    os.chdir(args.outdir) 
    numprefix=os.path.basename(bedfile_abpath)
    numfile=numprefix+".Numerical.txt"
    mapfile=numprefix+".map.txt"
    phenon=os.path.basename(pheno_abpath)

    pf=open(phenon,"r")
    line=pf.readline()
    pf.close()
    
    ls = line.split()
    trait=ls[1:]
    nfolder=len(trait)//farmcpujob +1
    lastfolder=len(trait)%farmcpujob

    


    k=2
    for i in range(1,nfolder):
        ed=k+(farmcpujob-1)
        os.system("mkdir xx.%d.%d-%d"%(i,k,ed))
        tmpr='''
        library(bigmemory)
        library(FarmCPUpp)
        setwd("%s")
        phnotype="%s"
        numgd="../%s"
        numif="../%s"

        myY <- read.table(phnotype,
                        header = TRUE, stringsAsFactors = FALSE)
        myGM <- read.table(numif,
                        header = TRUE, stringsAsFactors = FALSE)

        myGD <- read.big.matrix(numgd,
                                type = "double", sep = "\t", header = TRUE,
                                col.names = myGM$SNP, ignore.row.names = FALSE,
                                has.row.names = TRUE)
                                #, backingfile = "prunData_numeric.bin",
                                #descriptorfile = "prunData_numeric.desc")
        #dput(describe(myGD), "prunData_numeric_pointer.desc")

        # The big.matrix can be reattached in a different R session using
        # desc <- dget(pointer)
        # myGD <- attach.big.matrix(desc)
        for (i in %d:%d)
        {
        myResults <- farmcpu(Y = myY[,c(1,i)], GD = myGD, GM = myGM,ncores.glm=100 , ncores.reml=100)
        write_results(myResults)
        #manhattan_plot(myResults$COLEOPCOL_JI15$GWAS, cutoff = 0.01)
        rm(myResults)
        }

        '''%(os.getcwd()+"/"+"xx.%d.%d-%d"%(i,k,ed),pheno_abpath,numfile,mapfile,k,ed)

        outsh=open("xx.%d.%d-%d"%(i,k,ed)+"/"+"FarmCPU_%d_%d.r"%(k,ed),"w")
        outsh.write(tmpr)
        outsh.close()

        k+=farmcpujob
    
    if lastfolder >=1:
        os.system("mkdir xx.n.%d-%d"%(k,k+(lastfolder-1)))
        tmpr='''
        library(bigmemory)
        library(FarmCPUpp)
        setwd("%s")
        phnotype="%s"
        numgd="../%s"
        numif="../%s"

        myY <- read.table(phnotype,
                        header = TRUE, stringsAsFactors = FALSE)
        myGM <- read.table(numif,
                        header = TRUE, stringsAsFactors = FALSE)

        myGD <- read.big.matrix(numgd,
                                type = "double", sep = "\t", header = TRUE,
                                col.names = myGM$SNP, ignore.row.names = FALSE,
                                has.row.names = TRUE)
                                #, backingfile = "prunData_numeric.bin",
                                #descriptorfile = "prunData_numeric.desc")
        #dput(describe(myGD), "prunData_numeric_pointer.desc")

        # The big.matrix can be reattached in a different R session using
        # desc <- dget(pointer)
        # myGD <- attach.big.matrix(desc)
        for (i in %d:%d)
        {
        myResults <- farmcpu(Y = myY[,c(1,i)], GD = myGD, GM = myGM,ncores.glm=100 , ncores.reml=100)
        write_results(myResults)
        #manhattan_plot(myResults$COLEOPCOL_JI15$GWAS, cutoff = 0.01)
        rm(myResults)
        }

        '''%(os.getcwd()+"/"+"xx.n.%d-%d"%(k,k+(lastfolder-1)),pheno_abpath,numfile,mapfile,k,k+(lastfolder-1))

        outsh=open("xx.n.%d-%d"%(k,k+(lastfolder-1))+"/"+"FarmCPU_%d_%d.r"%(k,k+(lastfolder-1)),"w")
        outsh.write(tmpr)
        outsh.close()

        cmd='''
        for i in `ls | grep "xx\."`
        do
            %s $i/FarmCPU*.r &
        done
        wait 
        '''%(rpath )

        outsh=open("go.sh","w")
        outsh.write(cmd)
        outsh.close()

        #go farmcpu
        p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
        out,err = p.communicate()  #阻塞父进程
        for line in out.splitlines():
            print(line)
        for line in err.splitlines():
            print(line)

        #finish farmcpu
        cmd = "mkdir -p all_res;"
        os.system(cmd)
        os.chdir("./all_res")
        #ln data
        cmd='''
        ln -s ../xx*/*.csv .
        '''
        outsh=open("00.ln_data.sh","w")
        outsh.write(cmd)
        outsh.close()

        os.system(cmd)
        #plot
        plotr='''
        args=commandArgs(T)
        filen=args[1]
        titlenn=args[2]
        source("%s")
        #source("./GAPIT.QQ.FUNC_custom.r")
        source("%s")
        manhattanplot(mydata = filen, key = titlenn, chr = "all",title=titlenn, log10 = T, vline=F,x_tick_labs = "Character")
        mydata=read.table(filen,head=T)
        mydata=na.omit(mydata)
        jpeg(paste( titlenn,"_QQ-Plot.jpg" ,sep = ""),width = 15,height=15,units = 'cm', res = 225)
        qqunif.plot(mydata[,3],aspect=1)
        dev.off()
        #GAPIT.QQ(P.values = mydata[,3], name.of.trait = titlenn,DPP=50000)
        '''%(manplot,qqplot)

        outsh=open("plot.r","w")
        outsh.write(plotr)
        outsh.close()

        cmd='''
        for i in `ls *GWAS.csv`
        do
                ttn=`echo $i | awk -F ".GWAS" '{print $1}'`;
                echo "
                echo \\\"start...\\\"
                head -n 1 $i| awk  -F\\\",\\\"  'BEGIN{OFS=\\\"\\t\\\"}{print \$2,\$3,\$4}' > ${i}.addchr.csv;
                echo \\\"start generate plot data file...\\\"
                tail -n +2 $i | awk  -F\\\",\\\"  'BEGIN{OFS=\\\"\\t\\\"}{\$2=substr(\$1,1,5);print \$2,\$3,\$4}'  >> ${i}.addchr.csv;
                %s plot.r ${i}.addchr.csv ${ttn}
        " > $i.sh && chmod +x $i.sh && sh $i.sh &


        done
        wait
        '''%(rpath)

        outsh=open("01.plot.sh","w")
        outsh.write(cmd)
        outsh.close()

        os.system(cmd)

        #merge man plot and qq plot
        cmd = """
        #!/usr/bin/sh
        source /public/home/fengcong/anaconda2/etc/profile.d/conda.sh
        conda init
        conda activate py3
        %s %s %s


        conda deactivate
        """%(pypath,mergepy,"./")

        outsh=open("02.merge_jpg.sh","w")
        outsh.write(cmd)
        outsh.close()

        os.system("sh 02.merge_jpg.sh")





def go_gemma(bedfile,kinship,covf,maf,pheno,outdir,threads):
    # os.chdir(args.outdir) 
    #prepare genotype data and phenotype data
    bedfile_abpath=bedfile
    pheno_abpath=pheno
    kinship_abpath=kinship

    if bedfile.startswith("/"):
        #absolute path
        pass
    else:
        cwd=os.getcwd()
        bedfile_abpath = cwd+"/"+bedfile
    
    if pheno.startswith("/"):
        #absolute path
        pass
    else:
        cwd=os.getcwd()
        pheno_abpath = cwd+"/"+pheno
    if kinship:
        if kinship.startswith("/"):
            #absolute path
            pass
        else:
            cwd=os.getcwd()
            kinship_abpath = cwd+"/"+kinship

    cmd="""
    ln -s %s* %s;
    ln -s %s %s;
    """ %(bedfile_abpath,outdir,pheno_abpath,outdir)
    os.system(cmd)

    #threads
    if threads == "all":
        threads = 99999
    else:
        threads = int(threads)

    # generate script and run
    os.chdir(args.outdir) 
    bedn=os.path.basename(bedfile_abpath)
    phenon=os.path.basename(pheno_abpath)
    k = ""
    if kinship_abpath:
        k = "-k %s"%(kinship_abpath)
    else:
        k = "-k %s"%("./output/kinship.cXX.txt")
    
    
    cmd='''
    k=3
    n=1
    tb=0
    pho="%s"
    for i in `head -n 1 $pho  | cut -f 3- `
    do

            mkdir xx.${n}.${i}
            cd xx.${n}.${i}/;
            ln -s ../%s.bed .
            ln -s ../%s.bim .
            cut -f -5 ../%s.fam -d " "> tmp1.fam
            tail -n +2 ../$pho | cut -f $k  > tmp2.fam
            paste -d " " tmp1.fam tmp2.fam > %s.fam
            rm tmp*;
            cd -;

            let k+=1
            echo "
            %sgemma-0.98.1-linux-static -bfile %s -miss 0.9 -gk -o kinship 
            gemma-0.98.1-linux-static -bfile %s -miss 0.9 %s -lmm  %s %s -o GWAS
            " > xx.${n}.${i}/xx.${n}.${i}.gemma.sh && chmod +x xx.${n}.${i}/xx.${n}.${i}.gemma.sh  
            cd xx.${n}.${i}; 
            sh xx.${n}.${i}.gemma.sh 1>xx.${n}.${i}.gemma.sh.o 2>xx.${n}.${i}.gemma.sh.e & 
            cd -;
            let tb+=1
            if [ $tb -ge %d ];then
                let tb=0
                wait
            fi
            
            let n+=1

    done

    wait
    
    ''' %(phenon,bedn,bedn,bedn,bedn,"#" if kinship_abpath  else "",bedn,bedn,
        "-maf %s"%(maf) if maf else "" ,k, covf if covf=="" else " -c %s "%(covf) ,threads)

    outsh=open("0x.generate_script.sh","w")
    outsh.write(cmd)
    outsh.close()
    #exec 
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = p.communicate()  #阻塞父进程
    for line in out.splitlines():
        print(line)
    # p.wait()  

    #finish gemma 
    cmd = "mkdir -p all_res;"
    os.system(cmd)
    os.chdir("./all_res")

    cmd='''
    k=3
    n=1
    for i in `head -n 1 ../%s  | cut -f 3- `
    do
            #%s %s ../xx.${n}.${i}/output/GWAS.assoc.txt ../xx.${n}.${i}/*.fam ../xx.${n}.${i}/output/GWAS.assoc1.txt

            #ln -s ../xx.${n}.${i}/output/*assoc1.txt ./${i}.GWAS.assoc.txt
            ln -s ../xx.${n}.${i}/output/*assoc.txt ./${i}.GWAS.assoc.txt
            let k+=1
            let n+=1

    done

    wait
    '''%(phenon, rpath, padj)

    outsh=open("00.ln_data.sh","w")
    outsh.write(cmd)
    outsh.close()

    os.system(cmd)
    #plot 
    

    plotr='''
    args=commandArgs(T)
    filen=args[1]
    titlenn=args[2]
    source("%s")
    #source("./GAPIT.QQ.FUNC_custom.r")
    source("%s")
    manhattanplot(mydata = filen, key = titlenn, chr = "all",title=titlenn, log10 = T, vline=F,x_tick_labs = "Character")
    mydata=read.table(filen,head=T)
    mydata=na.omit(mydata)
    jpeg(paste( titlenn,"_QQ-Plot.jpg" ,sep = ""),width = 15,height=15,units = 'cm', res = 225)
    qqunif.plot(mydata[,3],aspect=1)
    dev.off()
    #GAPIT.QQ(P.values = mydata[,3], name.of.trait = titlenn,DPP=50000)
    '''%(manplot,qqplot)

    outsh=open("plot.r","w")
    outsh.write(plotr)
    outsh.close()



    cmd='''
    
    for i in `ls *GWAS.assoc.txt`
    do
            ttn=`echo $i | cut -f 1 -d "."`
            echo "
            echo \\\"start...\\\"
            head -n 1 $i| awk    'BEGIN{OFS=\\\"\\t\\\"}{print \$1,\$3,\$12}' > ${i}.addchr.csv;
            echo \\\"start generate plot data file...\\\"
            tail -n +2 $i | awk    'BEGIN{OFS=\\\"\\t\\\"}{print \$1,\$3,\$12}'  >> ${i}.addchr.csv;

            # echo \\\"start...\\\"
            # head -n 1 $i| awk    'BEGIN{OFS=\\\"\\t\\\"}{print \$1,\$3,\$13}' > ${i}.addchr.padj.csv;
            # echo \\\"start generate plot data file...\\\"
            # tail -n +2 $i | awk    'BEGIN{OFS=\\\"\\t\\\"}{\$1=substr(\$2,1,5);print \$1,\$3,\$13}'  >> ${i}.addchr.padj.csv;

            %s plot.r ${i}.addchr.csv ${ttn}
            # %s plot.r ${i}.addchr.padj.csv ${ttn}.padj
            
    " > $i.sh && chmod +x $i.sh && sh $i.sh &


    done
    wait
    ''' %(rpath,rpath)
    outsh=open("01.plot.sh","w")
    outsh.write(cmd)
    outsh.close()

    os.system(cmd)

    #merge man plot and qq plot
    cmd = """
    #!/usr/bin/sh
    source /public/home/fengcong/anaconda2/etc/profile.d/conda.sh
    conda init
    conda activate py3
    %s %s %s


    conda deactivate
    """%(pypath,mergepy,"./")

    outsh=open("02.merge_jpg.sh","w")
    outsh.write(cmd)
    outsh.close()

    os.system("sh 02.merge_jpg.sh")

    #get big than the threshold





 


def check_abs_path(filename):
    if filename.startswith("/"):
        return filename
    else:
        return os.getcwd()+"/" + filename

def write_sh(cmd,abs_filename):
    outsh=open(abs_filename,"w")
    outsh.write(cmd)
    outsh.close()

def go_emmax(tpedfile,phenotype,outdir,covf,threads):
    os.chdir(outdir)

    ##ln data
    cmd = '''
    ln -s %s.t* .
    ln -s %s .
    '''%(tpedfile,phenotype)
    write_sh(cmd,"00.ln_data.sh")
    os.system(cmd)

    cov=""
    if covf != "":
        os.system("ln -s %s ."%(covf))
        cov=os.path.basename(covf)

    #threads 
    if threads == "all":
        threads = 99999
    else:
        threads = int(threads)

    ##generate script and run
    tped=os.path.basename(tpedfile)
    phenon=os.path.basename(phenotype)
    
    cmd='''
    k=3
    n=1
    tb=0
    pho="%s"
    tpedprefix="%s"
    for i in `head -n 1 $pho  | cut -f 3- `
    do

            mkdir xx.${n}.${i}
            cd xx.${n}.${i}/;
            ln -s ../${tpedprefix}.tped .
            ln -s ../${tpedprefix}.tfam .
            tail -n +2 ../$pho|cut -f 1,2,$k  > ${tpedprefix}.pheno
            cd -;

            let k+=1
            echo "
            %s/emmax-kin-intel64 ${tpedprefix} -v -d 10 -o emmax_in.BN.kinf
            %s/emmax-intel64 -t ${tpedprefix} -o $i  -p ${tpedprefix}.pheno  -k emmax_in.BN.kinf
            " > xx.${n}.${i}/xx.${n}.${i}.emmax.sh && chmod +x xx.${n}.${i}/xx.${n}.${i}.emmax.sh  
            cd xx.${n}.${i}; 
            sh xx.${n}.${i}.emmax.sh 1>xx.${n}.${i}.emmax.sh.o 2>xx.${n}.${i}.emmax.sh.e & 
            #pid=$!
            #/public/agis/chengshifeng_group/fengcong/WGRS/software/cpulimit/cpulimit-0.2/bin/cpulimit -l 600 -p $pid -z -i &
            cd -;
            let tb+=1
            if [ $tb -ge %d ];then
                let tb=0
                wait
            fi
            let n+=1

    done

    wait
    
    ''' %(phenon,tped,emmaxpath,emmaxpath,threads)

    if covf != "":
        cmd='''
        k=3
        n=1
        tb=0
        pho="%s"
        tpedprefix="%s"
        for i in `head -n 1 $pho  | cut -f 3- `
        do

                mkdir xx.${n}.${i}
                cd xx.${n}.${i}/;
                ln -s ../${tpedprefix}.tped .
                ln -s ../${tpedprefix}.tfam .
                ln -s ../%s .
                tail -n +2 ../$pho|cut -f 1,2,$k  > ${tpedprefix}.pheno
                cd -;

                let k+=1
                echo "
                %s/emmax-kin-intel64 ${tpedprefix} -v -d 10 -o emmax_in.BN.kinf
                %s/emmax-intel64 -t ${tpedprefix} -o $i  -p ${tpedprefix}.pheno  -c %s -k emmax_in.BN.kinf
                " > xx.${n}.${i}/xx.${n}.${i}.emmax.sh && chmod +x xx.${n}.${i}/xx.${n}.${i}.emmax.sh  
                cd xx.${n}.${i}; 
                sh xx.${n}.${i}.emmax.sh 1>xx.${n}.${i}.emmax.sh.o 2>xx.${n}.${i}.emmax.sh.e & 
                cd -;
                let tb+=1
                if [ $tb -ge %d ];then
                    let tb=0
                    wait
                fi
                let n+=1

        done

        wait
        
        ''' %(phenon,tped,cov,emmaxpath,emmaxpath,cov,threads)

    write_sh(cmd,"01.generate_script.sh")
    os.system(cmd)

    ## finish emmax 
    cmd = "mkdir -p all_res;"
    os.system(cmd)
    os.chdir("./all_res")

    cmd='''
    k=3
    n=1
    for i in `head -n 1 ../%s  | cut -f 3- `
    do


            ln -s ../xx.${n}.${i}/*.ps ./
            let k+=1
            let n+=1

    done
    '''%(phenon)

    write_sh(cmd,"00.ln_data.sh")
    os.system(cmd)

    #
    plotr='''
    args=commandArgs(T)
    filen=args[1]
    titlenn=args[2]
    source("%s")
    #source("./GAPIT.QQ.FUNC_custom.r")
    source("%s")
    manhattanplot(mydata = filen, key = titlenn, chr = "all",title=titlenn, log10 = T, vline=F,x_tick_labs = "Character")
    mydata=read.table(filen,head=T)
    mydata=na.omit(mydata)
    jpeg(paste( titlenn,"_QQ-Plot.jpg" ,sep = ""),width = 15,height=15,units = 'cm', res = 225)
    qqunif.plot(mydata[,3],aspect=1)
    dev.off()
    #GAPIT.QQ(P.values = mydata[,3], name.of.trait = titlenn,DPP=50000)
    '''%(manplot,qqplot)

    write_sh(plotr,"plot.r")

    #prepare plot data
    cmd='''
    
    for i in `ls *.ps`
    do
            ttn=`echo $i | cut -f 1 -d "."`
            echo "
            echo \\\"start...\\\"
            echo -ne \\\"chr\\tpos\\tpvalue\\n\\\" > ${i}.addchr.csv;
            echo \\\"start generate plot data file...\\\"
            awk    'BEGIN{OFS=\\\"\\t\\\"}{\$2=substr(\$1,1,5);\$3=substr(\$1,7);print \$2,\$3,\$4}'  $i >> ${i}.addchr.csv;
            #cp plot.r ${i}.addchr.csv.r && sed -i \\\"s:filen:${i}.addchr.csv:g\\\" ${i}.addchr.csv.r && sed -i \\\"s:titlenn:${ttn}:g\\\" ${i}.addchr.csv.r
            #/gxyy-local/home/liyongyao/fengcong/software/R/bin/R CMD BATCH --args ${i}.addchr.csv.r 
            %s plot.r ${i}.addchr.csv ${ttn}
            
    " > $i.sh && chmod +x $i.sh && sh $i.sh &


    done
    wait
    ''' %(rpath)

    write_sh(cmd,"01.plot.sh")
    os.system(cmd)

    #merge man plot and qq plot
    cmd = """
    #!/usr/bin/sh
    source /public/home/fengcong/anaconda2/etc/profile.d/conda.sh
    conda init
    conda activate py3
    %s %s %s


    conda deactivate
    """%(pypath,mergepy,"./")

    write_sh(cmd,"02.merge_jpg.sh")
    os.system("sh 02.merge_jpg.sh")

def go_cmlm_mlm(hmpfile,phenotype,outdir,threads,mod):
    os.chdir(outdir)

    ##ln data
    cmd = '''
    ln -s %s .
    ln -s %s .
    '''%(hmpfile,phenotype)
    write_sh(cmd,"00.ln_data.sh")
    os.system(cmd)

    #threads 
    if threads == "all":
        threads = 99999
    else:
        threads = int(threads)

    ##generate script and run
    hmp=os.path.basename(hmpfile)
    phenon=os.path.basename(phenotype)

    cmd='''
    k=3
    n=1
    tb=0
    pho="%s"
    hmp="%s"
    mod=%s
    for i in `head -n 1 $pho  | cut -f 3- `
    do

            mkdir xx.${n}.${i}
            cd xx.${n}.${i}/;
            ln -s ../${hmp} .
            tail -n +2 ../$pho|cut -f 2,$k  > ${hmp}.pheno
            cd -;

            let k+=1
            echo '

            library(multtest)
            library(gplots)
            library(LDheatmap)
            library(genetics)
            library(ape)
            library(EMMREML)
            library(compiler) #this library is already installed in R library("scatterplot3d")
            library("scatterplot3d")
            source("/public/agis/chengshifeng_group/fengcong/WGRS/software/gwas_pipline/gapit_functions.txt")
            source("/public/agis/chengshifeng_group/fengcong/WGRS/software/gwas_pipline/emma.txt")
            myY <- read.table("%s.pheno", head = F)
            myG <- read.table("%s", head = FALSE,comment.char = "")

            #Step 2: Run GAPIT (CMLM Model)
            myGAPIT <- GAPIT(
            Y=myY,
            G=myG,
            Model.selection = TRUE,
            model="%s",
            PCA.total=3
            )


            ' > xx.${n}.${i}/xx.${n}.${i}.${mod}.R  
            cd xx.${n}.${i}; 

            %s xx.${n}.${i}.${mod}.R 1>xx.${n}.${i}.${mod}.R.o 2>xx.${n}.${i}.${mod}.R.e & 
            cd -;
            let tb+=1
            if [ $tb -ge %d ];then
                let tb=0
                wait
            fi
            let n+=1

    done

    wait
    
    ''' %(phenon,hmp,mod,hmp,hmp,mod,rpath,threads)

    write_sh(cmd,"01.generate_script.sh")
    os.system(cmd)

    ## finish cmlm 
    cmd = "mkdir -p all_res;"
    os.system(cmd)
    os.chdir("./all_res")

    #ln data
    cmd='''
    k=3
    n=1
    for i in `head -n 1 ../%s  | cut -f 3- `
    do


            ln -s ../xx.${n}.${i}/*.GWAS.Results.csv ./$i.GWAS.Results.csv
            let k+=1
            let n+=1

    done
    ''' %(phenon)
    write_sh(cmd,"00.ln_data.sh")
    os.system(cmd)

    #
    plotr='''
    args=commandArgs(T)
    filen=args[1]
    titlenn=args[2]
    source("%s")
    #source("./GAPIT.QQ.FUNC_custom.r")
    source("%s")
    manhattanplot(mydata = filen, key = titlenn, chr = "all",title=titlenn, log10 = T, vline=F,x_tick_labs = "Character")
    mydata=read.table(filen,head=T)
    mydata=na.omit(mydata)
    jpeg(paste( titlenn,"_QQ-Plot.jpg" ,sep = ""),width = 15,height=15,units = 'cm', res = 225)
    qqunif.plot(mydata[,3],aspect=1)
    dev.off()
    #GAPIT.QQ(P.values = mydata[,3], name.of.trait = titlenn,DPP=50000)
    '''%(manplot,qqplot)

    write_sh(plotr,"plot.r")

    #prepare plot data
    cmd='''
    
    for i in `ls *.csv`
    do
            ttn=`echo $i | cut -f 1 -d "."`
            echo "
            echo \\\"start...\\\"
            echo -ne \\\"chr\\tpos\\tpvalue\\n\\\" > ${i}.addchr.csv;
            echo \\\"start generate plot data file...\\\"
            tail -n +2 $i| awk  -F \\\",\\\"  'BEGIN{OFS=\\\"\\t\\\"}{\$2=substr(\$1,1,5);\$3=substr(\$1,7);print \$2,\$3,\$4}'   >> ${i}.addchr.csv;
            %s plot.r ${i}.addchr.csv ${ttn}
            
    " > $i.sh && chmod +x $i.sh && sh $i.sh &


    done
    wait
    ''' %(rpath)

    write_sh(cmd,"01.plot.sh")
    os.system(cmd)

    cmd = """
    #!/usr/bin/sh
    source /public/home/fengcong/anaconda2/etc/profile.d/conda.sh
    conda init
    conda activate py3
    %s %s %s


    conda deactivate
    """%(pypath,mergepy,"./")

    write_sh(cmd,"02.merge_jpg.sh")
    os.system("sh 02.merge_jpg.sh")


if __name__ == "__main__":
    START_TIME = datetime.now()

    cmdparser = argparse.ArgumentParser(description="GWAS pipline(gemma/farmcpupp/emmax)")
    cmdparser.add_argument("-t", "--tped", dest="tped", type=str, 
                           help="prefix of plink .tped file .")
    cmdparser.add_argument("-c", "--covf", dest="covf", type=str, 
                           help="multi-column covariate file with FAMID, INDID at the first two colmns . if model == gamma ,the first colum must be 1, and dont need FID and IID.")
    cmdparser.add_argument("-b", "--bed", dest="prefix", type=str, 
                           help="prefix of plink .bed file .")
    cmdparser.add_argument("-k", "--kinship", dest="kinship", type=str, 
                           help="kinship file for gemma.")
    cmdparser.add_argument("-f", "--maf", dest="maf",type=str, 
                           help="maf for gemma.")
    # cmdparser.add_argument("-v","--vcf", dest="vcf",type=str, 
    #                         help="input vcf file.")
    cmdparser.add_argument("-n","--numprefix", dest="numprefix",type=str, 
                             help="input numprefix.Numerical.txt file,used for FarmCPUpp.")
    cmdparser.add_argument("-m","--hmp", dest="hmp",type=str, 
                             help="hmp file,used for MLM and CMLM.50GB mem per 1Mb marker.")
    cmdparser.add_argument("-p","--pheno", dest="pheno",type=str, required=True,
                            help="phenotype data.")
    cmdparser.add_argument("-o","--outdir", dest="outdir",type=str, required=True,
                            help="output dir.")
    cmdparser.add_argument("-@","--threads", dest="threads",type=str,default="all",
                            help="How many traits do we compute in parallel.default:all")
    cmdparser.add_argument("-d","--model", dest="model",type=str, required=True,
                            help="model choose:gemma/farmcpu/emmax/mlm/cmlm")
    
    args = cmdparser.parse_args()

    #args.bed
    print(args)
    ##mkdir 
    cmd='''
    mkdir -p %s;
    ''' %(args.outdir)
    os.system(cmd)

    ##choose model
    md={"gemma","farmcpu","emmax","cmlm","mlm"}
        
    if args.model.lower() == "gemma":
        bedfile=args.prefix
        phenotype=args.pheno
        kinship = args.kinship
        covf=""
        if args.covf != None:
            covf=check_abs_path(args.covf)
        maf=args.maf
        
        go_gemma(bedfile,kinship,covf,maf,phenotype,args.outdir,args.threads)
    elif args.model.lower() == "farmcpu":
        numfile=args.numprefix
        phenotype=args.pheno
        go_farmcpupp(numfile,phenotype,args.outdir)
    elif args.model.lower() == "emmax":
        tpedfile=check_abs_path(args.tped)
        phenotype=check_abs_path(args.pheno)
        od=check_abs_path(args.outdir)
        covf=""
        if args.covf != None:
            covf=check_abs_path(args.covf)

        go_emmax(tpedfile,phenotype,od,covf,args.threads)
    elif args.model.lower() == "cmlm" or args.model.lower() == "mlm":
        hmpfile=check_abs_path(args.hmp)
        phenotype=check_abs_path(args.pheno)
        od=check_abs_path(args.outdir)
        go_cmlm_mlm(hmpfile,phenotype,od,args.threads,args.model.upper())
    else:
        print("model choose error:",md)
        exit(-1)
