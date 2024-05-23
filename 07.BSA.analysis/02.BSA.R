library("vcfR")

vcf <- read.vcfR("/vol3/agis/chengshifeng_group/shiyan/06-BSA_seq/03.fa/fa.SNP.HARD.Missing-unphasing.ID.allele2_retain.hard_retain.vcf")

chrom <- getCHROM(vcf)
pos <- getPOS(vcf)
ref <- getREF(vcf)
alt <- getALT(vcf)

ad <- extract.gt(vcf, "AD")
ref_split <- masplit(ad, record = 1, sort = 0)
alt_split <- masplit(ad, record = 2, sort = 0)
gt <- extract.gt(vcf, "GT")

df <- data.frame(CHROM = chrom,
                 POS = pos,
                 REF = ref,
                 ALT = alt,
                 AD_REF.GP_G = ref_split[,1],
                 AD_ALT.GP_G = alt_split[,1],
                 AD_REF.GP_Y = ref_split[,2],
                 AD_ALT.GP_Y = alt_split[,2]
)

mask <- which(gt[,"JI2822"] != "0/1" &  gt[,"JI0816"] != "0/1")

df1 <- df[mask,]


write.table(df1, file = "GP1.tsv", sep = "\t", row.names = F, quote = F)

#######BSA分析以及可视化
setwd("D:/ChengLab/豌豆/七大性状/七大性状整理/花的位置/")
library(QTLseqr)
library(ggplot2)

df<-importFromTable("Fa.tsv",
                    highBulk = 'Fa_fa',
                    lowBulk = 'Fa_WT',
                    chromList = paste0('chr',formatC(1:7,width = 1,flag =0)),
                    sep = "\t")
df<-subset(df,!is.na(SNPindex.LOW)&!is.na(SNPindex.HIGH))
#通常还会增加一个SNP的过滤，主要是结合数据具体深度和等位基因频率来进行

#深度，通过查看分布，设置临界值为80
depth_plot<-ggplot(data=df)+
  geom_histogram(aes(x=DP.HIGH+DP.LOW))+
                     xlim(0,200)

# 等位基因频率，0.2-0.8
alle_plot <- ggplot(data = df) +           
  geom_histogram(aes(x = REF_FRQ))

#以高值池SNP-index分布，低值池类似，通常情况，大部分的SNP应该符合正态分布，在0.5左右，两侧的极端通常是在两个混池中基因型一致的样本，可以剔除。通常根据这个检查数据是否异常
high_index_plot <- ggplot(data = df) +
  geom_histogram(aes(x = SNPindex.LOW))

df <-filterSNPs(
  SNPset = df,
  refAlleleFreq = 0.20, #等位基因频率过滤，0.2表示0.2~0.8之间
  minTotalDepth = 10, #最小过滤深度
  maxTotalDepth = 120, #最大过滤深度
#  minSampleDepth = 40, #单个样本过滤深度
#  minGQ = 99, ### genotype quality 过滤
  verbose = TRUE
)



#G statis
df<-runGprimeAnalysis(SNPset=df,windowSize = 1e+07,outlierFilter = "deltaSNP",filterThreshold = 0.1)

#delta SNP 置信区间
df<-runQTLseqAnalysis(SNPset = df,windowSize = 1e7,popStruc="F2",bulkSize = c(20,20))


write.csv(df,"Fa.BSA.tsv",sep="\t",row.names = T,col.names = T)

#plot

g_plot<-plotQTLStats(SNPset = df,var="Gprime",plotThreshold = TRUE,q=0.01)+theme_bw()
delta_plot<-plotQTLStats(SNPset = df, var="deltaSNP",plotIntervals = TRUE)+theme_bw()

pdf("Fa.camor.BSA.plot.pdf",height = 8,width = 21)
g_plot
delta_plot
dev.off()
