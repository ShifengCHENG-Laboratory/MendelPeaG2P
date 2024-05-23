manhattanplot <- function(
  mydata,					#输入文件
  key,						#输出文件的前缀
  chr = "all",				#如果选择染色体n，则为chr = "n"
  start = NA,				# 绘图的起始位置
  end = NA,
  type = "point",			#可以多选，type = c("point", "bar", "line")
  base_value = "min",		# 绘图的基线，Y轴如果有负值，则用默认值base_value = "min"。Y轴值有小于0，则设置为base_value = 0 。
  reflen = NULL,
  columns=c(1,2,3),			#选择绘图数据所在的列，分别为染色体iD，位置和Y轴数值。
  log10 = TRUE,				#是否数值取log10
  log2 = FALSE,
  zscore = FALSE,			#是否Z-score处理Y轴数值
  minus = FALSE,
  vline = TRUE,				#染色体之间是否加垂直虚线
  output_plot = TRUE,
  maincol = c("#808080","#000000"),            #c("#db1b2b", "#72c334", "#db1b2b", "#413496", "#495226", "#d60b6f", "#e66519", "#4197d8", "#d581b7", "#83d3ad", "#7c162c", "#26755d"),
  cexMain = 15,
  cexTick = 15,
  title = NULL,
  threshold_line_anyway = TRUE,
  threshold_line_size = 0.6,
  thresholds = c(-log10(0.01/nrow(na.omit(mydf))), -log10(0.1/nrow(na.omit(mydf)))),		 	#阈值线，默认两条
  thresholds_color = c("red", "blue","red"),															#阈值线的颜色，对应默认定两个色
  x_tick = TRUE,
  x_tick_labs = "Numeric",               # 默认x label只写数字。如果染色体中存在字符且期望显示字符，则改为= "character"                
  x_tick_angle = 0,
  x_tick_vjust = 0.5,
  xlab = "Chromosomes",
  ylab = expression(-log[10](italic(p))),				# Y轴标签
  ymin = NULL,
  ymax = NULL,
  pointSize = 0.8,
  width = 9*2.5,
  height = 4.5*2.5
){
  ############### function for multipple mixed order ###############
  ### loading packages
  library(gtools)
  library(ggplot2)
  library(scales)
  ### multiple mixed order function
  multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
    do.call(order, c(
      lapply(list(...), function(l){
        if(is.character(l)){
          factor(l, levels=mixedsort(unique(l)))
        } else {
          l
        }
      }),
      list(na.last = na.last, decreasing = decreasing)
    ))
  }
  ### read data
  mydf <- read.table(mydata, header = T, stringsAsFactors = F, comment.char = "", check.names = F)
  mydf <- mydf[, columns]
  mydf[,2] <- as.numeric(mydf[,2])
  print("read in columns:")
  print(head(mydf))
  names(mydf) <- c("Chr", "Pos", "P")
  mydf$Chr <- as.character(mydf$Chr)
  chr <- as.character(chr)
  ### log trans
  if(log10 == T) {
    mydf$P <- -log10(mydf$P)
  }
  if(log2 == T) {
    mydf$P <- -log2(mydf$P)
  }
  if(zscore == T) {
    mydf$P <- scale(mydf$P)
  }
  if(minus == T) {
    mydf$P <- mydf$P * -1
  }
  ### threshold value
  mydf <- mydf[!is.infinite(mydf[, 3]), ]
  mydf <- na.omit(mydf)
  mydf <<- mydf
  gcl_and_sgl <- vector()
  for (i in 1:length(thresholds)) {
    gcl_and_sgl[i] <- as.numeric(thresholds[i])
  }
  gcl_and_sgl <<- gcl_and_sgl
  gcl <- max(gcl_and_sgl)
  if(length(gcl_and_sgl) != length(thresholds_color)){
    thresholds_color <- rep(thresholds_color, length(gcl_and_sgl))[1:length(gcl_and_sgl)]
  }
  print("threshold:")
  print(gcl_and_sgl)
  print(thresholds_color)
  chrID <- mixedsort(unique(mydf[,1]))
  ###
  if(!is.null(reflen)){
    myreflen <- read.table(reflen, header = F, stringsAsFactors = F, comment.char = "")
    chr <- NULL
  }else if(!is.null(chr) & is.null(reflen)){
    if ("all" %in% chr){
      myreflen <- data.frame(chr = chrID, stringsAsFactors = F)
    } else {
      myreflen <- data.frame(chr = chr, stringsAsFactors = F)
    }
  } else if(!is.null(chr) & !is.null(reflen)){
    print("Use reflen when reflen and chr are all defined!\n")
    myreflen <- read.table(reflen, header = F, stringsAsFactors = F, comment.char = "")
    chr = NULL
  }
  myreflen <- data.frame(myreflen[multi.mixedorder(myreflen[,1]),], stringsAsFactors = F)
  mydf <- mydf[mydf[,1] %in% myreflen[,1],]
  if(!is.na(start) & !is.na(end)){
    mydf <- subset(mydf, Pos >= start & Pos <= end)
  }
  print("sub:")
  print(head(mydf))
  ############## deal with chromosome IDs and make sure positions in a series #############
  chrDF <- data.frame(OriginID = myreflen[, 1], NumericID = 1:nrow(myreflen), stringsAsFactors = F)
  chrNum <- nrow(chrDF)
  write.table(chrDF, file = paste0(key, "_chrID_table.txt"), row.names = F, quote = F, col.names = T, sep = "\t")
  mydf <- mydf[multi.mixedorder(mydf$Chr,mydf$Pos),]
  print("ordered:")
  print(head(mydf))
  print(chrNum)
  for(i in 1:chrNum) {
    mydf$New_Chr[mydf$Chr == chrDF[i,1]] <- chrDF[i,2]
    index <- which(mydf[, 1] == chrDF[i, 1])
    lastmarker <- max(mydf[index, 2])
    ##print(lastmarker)
    if (i < chrNum){
      index2 <- which(mydf[, 1] == chrDF[i + 1, 1])
      mydf[index2, 2] <- mydf[c(index2), 2] + lastmarker
    }
  }
  ############## middle and max of each chromosome ############
  bpMid <- vector(length = chrNum)
  bpMax <- vector(length = chrNum)
  if(chrNum == 1){
    if(max(mydf$Pos) - min(mydf$Pos) > 2000000){
      mydf$Pos <- mydf$Pos/1000000
      xlab = "Positions(Mb)"
    }else if(max(mydf$Pos) - min(mydf$Pos) > 5000 & max(mydf$Pos) - min(mydf$Pos) <= 2000000){
      mydf$Pos <- mydf$Pos/1000
      xlab = "Positions(Kb)"
    }
  }
  for(i in 1:chrNum){
    index <- which(mydf[, 1] == chrDF[i, 1])
    posSub <- mydf[index, 2]
    bpMax[i] <- max(posSub)
    bpMid[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
  }
  ############### create manhattan plot ##############
  if(maincol[1] == "rainbow") {
    maincol <- rainbow(chrNum)
  } else {
    maincol <- rep(maincol, chrNum)[1:chrNum]
  }
  print("max P value:")
  print(max(mydf$P))
  if(is.null(ymax)){
    ymax = max(mydf$P)*1.1
  }else{
    if(max(mydf$P) > ymax){
      ymax = ymax
      mydf$P[mydf$P > ymax] <- ymax - 2;
    }else{
      ymax = max(mydf$P)
    }
  }
  if(is.null(ymin)){
    ymin = min(mydf$P)
  }else{
    ymin = ymin
  }
  if(ymax <= gcl) {
    ymax = gcl*1.1
  }
  print(head(mydf))
  ############################################
  if("point" %in% type){
    ##### original plot #####
    if(vline == T & chrNum > 1) {
      p <- ggplot(mydf) + geom_vline(xintercept = bpMax, linetype = 1, col = "grey")
    } else {
      p <- ggplot(mydf)
    }
    p <- p + geom_point(aes(x = Pos, y = round(P, digits = 20), colour = as.factor(New_Chr), size = P), shape = 20)
    p <- p + scale_size_continuous(range = c(pointSize, pointSize*1.8))
    ##### decide x tick label type(numeric or characters) #####
    if(x_tick == T) {
      ##################### with x ticks ########################
      if (x_tick_labs == "Numeric") {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 2]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- 0
        x_tick_vjust <- 0
      } else {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 1]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- x_tick_angle
        x_tick_vjust <- x_tick_vjust
      }
    } else {
      ############### without x ticks ####################
      p <- p + scale_x_continuous(breaks = NULL, expand = c(0.005, 0))
    }
    p <- p + scale_y_continuous(expand = c(0.01,0), limits = c(ymin, ymax), breaks = pretty_breaks())
    #set themes
    p <- p + scale_color_manual(values = c(maincol)) + theme_classic()
    ##### xlab and ylab #####
    if(!is.null(xlab)){
      p <- p + xlab(xlab)
    } else {
      p <- p + theme(axis.title.x=element_blank())
    }
    if(!is.null(ylab)){
      p <- p + ylab(ylab)
    } else {
      p <- p + theme(axis.title.y=element_blank())
    }
    ##### threshold lines ####
    if(sum(is.na(gcl_and_sgl)) == 0){
      line_df <- data.frame(threshold = c(line = gcl_and_sgl))
      p <- p + geom_hline(data = line_df, aes(yintercept = threshold), size = threshold_line_size, linetype = 2, col = thresholds_color)
    }
    ##### final themes #####
    p <- p + theme(legend.position = "none",
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.border = element_blank(),
                   panel.background = element_blank(),
                   axis.title.x = element_text(face = "bold", size = cexMain, colour = "black"),
                   axis.title.y = element_text(face = "bold", size = cexMain, colour = "black"),
                   axis.text.x = element_text(size = cexTick,angle = x_tick_angle, vjust = x_tick_vjust, face = "bold", colour = "black"), 
                   axis.text.y = element_text(face = "bold", size = cexTick, colour = "black"),
                   axis.line.x = element_line(colour = "black", size = 1),
                   axis.line.y = element_line(colour = "black", size = 1))
    if(!is.null(title)){
      p <- p + ggtitle(title)
      p <- p + theme(plot.title = element_text(size=cexMain, face="bold", hjust = 0.5))
    }
    #return(p)
    if(output_plot == T) {
      #pdf(file = paste0(key, "_manhattan.pdf"), width = width/2.45, height = height/2.45)
      #print(p)
      #dev.off()
      # png(paste0(key,"_manhattan.png"), width = width, height = height, units = 'cm', res = 300)
      # print(p)
      # dev.off()
      jpeg(paste0(key,"_manhattan.jpg"), width = width, height = height, units = 'cm', res = 300)
      print(p)
      dev.off()
      #    svg(filename = paste0(key, "_manhattan.svg"), width = width/2.54, height = height/2.54)
      #    print(p)
      #    dev.off()
    }
  }
  if("bar" %in% type){
    #mydf$New_Chr <- factor(mydf$New_Chr, levels = order(as.numeric(unique(mydf$New_Chr))))
    p <- ggplot(mydf)
    if(vline == T){
      p <- p + geom_vline(xintercept = bpMax, linetype = 1, col = "grey")
    }
    if(base_value == "min"){
      p <- p + geom_rect(aes(xmin = Pos, xmax = Pos,
                             ymin = min(P), ymax = P,
                             fill= as.factor(New_Chr), colour = as.factor(New_Chr)),
                         position = "identity", stat = "identity")
    } else if(is.numeric(base_value)){
      p <- p + geom_rect(aes(xmin = Pos, xmax = Pos,
                             ymin = base_value, ymax = P,
                             fill= as.factor(New_Chr), colour = as.factor(New_Chr)),
                         position = "identity", stat = "identity")
    } else {
      print("base_value must be min or numeric")
      print(base_value)
      print("is not min or any numeric value")
    }

    if(x_tick == T) {
      ##################### with x ticks ########################
      if (x_tick_labs == "Numeric") {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 2]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- 0
        x_tick_vjust <- 0
      } else {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 1]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- x_tick_angle
        x_tick_vjust <- x_tick_vjust
      }
    } else {
      ############### without x ticks ####################
      p <- p + scale_x_continuous(breaks = NULL, expand = c(0.005, 0))
    }
    p <- p + scale_y_continuous(expand = c(0.01,0), limits = c(ymin, ymax), breaks = pretty_breaks())
    #set themes
    p <- p + scale_color_manual(values = c(maincol)) + theme_classic()
    ##### xlab and ylab #####
    if(!is.null(xlab)){
      p <- p + xlab(xlab)
    } else {
      p <- p + theme(axis.title.x=element_blank())
    }
    if(!is.null(ylab)){
      p <- p + ylab(ylab)
    } else {
      p <- p + theme(axis.title.y=element_blank())
    }
    p <- p + theme_classic() + theme(legend.position = "none",
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.border = element_blank(),
                                     panel.background = element_blank(),
                                     axis.title.x = element_text(face = "bold", size = cexMain, colour = "black"),
                                     axis.title.y = element_text(face = "bold", size = cexMain, colour = "black"),
                                     axis.text.x = element_text(size = cexTick,angle = x_tick_angle, vjust = x_tick_vjust, face = "bold", colour = "black"), 
                                     axis.text.y = element_text(face = "bold", size = cexTick, colour = "black"),
                                     axis.line.x = element_line(colour = "black", size = 1),
                                     axis.line.y = element_line(colour = "black", size = 1))
    if(sum(is.na(gcl_and_sgl)) == 0){
      line_df <- data.frame(threshold = c(line = gcl_and_sgl))
      p <- p + geom_hline(data = line_df, aes(yintercept = threshold), size = threshold_line_size, linetype = 2, col = thresholds_color)
    }
    if(!is.null(title)){
      p <- p + ggtitle(title)
      p <- p + theme(plot.title = element_text(size=cexMain, face="bold", hjust = 0.5))
    }
    #pdf(file = paste0(key, "_barplot.pdf"), width = width/2.54, height = height/2.54)
    #print(p)
    #dev.off()
    # png(paste0(key,"_barplot.png"), width = width, height = height, units = 'cm', res = 300)
    # print(p)
    # dev.off()
    jpeg(paste0(key,"_barplot.jpg"), width = width, height = height, units = 'cm', res = 300)
    print(p)
    dev.off()
  }
  if("line" %in% type){
    p <- ggplot(mydf)
    if(vline == T){
      p <- p + geom_vline(xintercept = bpMax, linetype = 1, col = "grey")
    }
    p <- p + geom_line(aes(x = Pos, y = P, color = as.factor(New_Chr)))
    if(x_tick == T) {
      ##################### with x ticks ########################
      if (x_tick_labs == "Numeric") {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 2]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- 0
        x_tick_vjust <- 0
      } else {
        if(length(bpMid) == 1){
          ##### one chromosome only #######
          p <- p + scale_x_continuous(breaks = pretty_breaks(n = 5), expand = c(0.005, 0))
        }else{
          p <- p + scale_x_continuous(labels = as.character(chrDF[, 1]), breaks = bpMid, expand = c(0.005, 0))
        }
        x_tick_angle <- x_tick_angle
        x_tick_vjust <- x_tick_vjust
      }
    } else {
      ############### without x ticks ####################
      p <- p + scale_x_continuous(breaks = NULL, expand = c(0.005, 0))
    }
    p <- p + scale_y_continuous(expand = c(0.01,0), limits = c(ymin, ymax), breaks = pretty_breaks())
    #set themes
    p <- p + scale_color_manual(values = c(maincol)) + theme_classic()
    ##### xlab and ylab #####
    if(!is.null(xlab)){
      p <- p + xlab(xlab)
    } else {
      p <- p + theme(axis.title.x=element_blank())
    }
    if(!is.null(ylab)){
      p <- p + ylab(ylab)
    } else {
      p <- p + theme(axis.title.y=element_blank())
    }
    p <- p + theme_classic() + theme(legend.position = "none",
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.border = element_blank(),
                                     panel.background = element_blank(),
                                     axis.title.x = element_text(face = "bold", size = cexMain, colour = "black"),
                                     axis.title.y = element_text(face = "bold", size = cexMain, colour = "black"),
                                     axis.text.x = element_text(size = cexTick,angle = x_tick_angle, vjust = x_tick_vjust, face = "bold", colour = "black"), 
                                     axis.text.y = element_text(face = "bold", size = cexTick, colour = "black"),
                                     axis.line.x = element_line(colour = "black", size = 1),
                                     axis.line.y = element_line(colour = "black", size = 1))
    if(sum(is.na(gcl_and_sgl)) == 0){
      line_df <- data.frame(threshold = c(line = gcl_and_sgl))
      p <- p + geom_hline(data = line_df, aes(yintercept = threshold), size = threshold_line_size, linetype = 2, col = thresholds_color)
    }
    if(!is.null(title)){
      p <- p + ggtitle(title)
      p <- p + theme(plot.title = element_text(size=cexMain, face="bold", hjust = 0.5))
    }
    #pdf(file = paste0(key, "_lineplot.pdf"), width = width/2.54, height = height/2.54)
    #print(p)
    #dev.off()
    # png(paste0(key,"_lineplot.png"), width = width, height = height, units = 'cm', res = 300)
    # print(p)
    # dev.off()
    jpeg(paste0(key,"_lineplot.jpg"), width = width, height = height, units = 'cm', res = 300)
    print(p)
    dev.off()
  }
}
