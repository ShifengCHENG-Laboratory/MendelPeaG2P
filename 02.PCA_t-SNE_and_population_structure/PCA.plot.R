pca_plot <- function(
  eigenvector,
  group,
  outdir,
  key = "pca",
  size = 3,
  shape = F,
  shapes = NULL,
  width = 8.5,
  height = 6,
  cexMain = 18,
  cexTick = 15,
  x_tick_angle = 0,
  x_tick_vjust = 0,
  add_lab = F,
  lab_sample_list = c(),
  expand = c(0.02, 0.02),
  border = F,
  border_size = 2,
  line0 = F,
  line0_size = 0.5,
  line0_color = "grey",
  alpha = 1
){
  library(ggplot2)
  library(ggrepel)
  library(scales)
  setwd(outdir)
  the <- theme(#legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = cexMain, face = "bold.italic"),
    legend.key.height=unit(2,"line"),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 1, colour = "black"),
    axis.ticks.length=unit(.15, "cm"),
    axis.title.x = element_text(face = "bold", size = cexMain, colour = "black"),
    axis.title.y = element_text(face = "bold", size = cexMain, colour = "black"),
    axis.text.x = element_text(size = cexTick,angle = x_tick_angle, vjust = x_tick_vjust, face = "bold", colour = "black"), 
    axis.text.y = element_text(face = "bold", size = cexTick, colour = "black"))
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1))
  egv <- read.table(eigenvector, header = T, stringsAsFactors = F,sep = "\t")
  #eva <- read.table(eigenvalue, header = T, stringsAsFactors = F)
  pop <- read.table(group, header = T, stringsAsFactors = F, check.names = F, comment.char = "", sep = "\t")
  egv_order <- egv[match(pop[,3], egv[,1]),]
  pop_order <- pop[match(egv_order[,1], pop[,3]),]
  egv_order$group <- pop_order[,4]
  egv_order$group <- as.factor(egv_order$group)
  if(ncol(pop_order) >= 5){
    egv_order$color <- pop_order[,5]
    col_df <- unique(pop_order[,4:5])
    col_df <- col_df[order(col_df[,1]),]
    print(col_df)
  }
  ### main plot
  p12 <- ggplot(data = egv_order, aes(x = Dim1, y = Dim2, group = group, color = group))
  
  
  if(line0 == T){
    p12 <- p12 + geom_hline(yintercept = 0, linetype = "dashed", color = line0_color, size = line0_size)
    p12 <- p12 + geom_vline(xintercept = 0, linetype = "dashed", color = line0_color, size = line0_size)
    
  }
  
  if(shape == T){
    p12 <- p12 + geom_point(aes(shape = group, fill = group), size = size, alpha = alpha)

  }else{
    p12 <- p12 + geom_point(size = size, alpha = alpha, shape = 19)

  }
  
  if(add_lab == T){
    key <- paste0(key, ".label")
    egv_order$label <- ""
    egv_order[match(lab_sample_list, egv_order[,1]),]$label <- lab_sample_list
    p12 <- p12 + geom_text_repel(data = egv_order , vjust = -1, aes(x = Dim1, y = Dim2, label = label),show.legend = F,box.padding = 2, max.overlaps =600)

  }
  
  if(ncol(pop_order) >= 5){
    p12 <- p12 + scale_color_manual(labels = col_df$group, values = col_df$color)
    p12 <- p12 + scale_fill_manual(values = col_df$color)
    

  }
  
  ### PC1 and PC2 plot
  if(!is.null(shapes)){
    p12 <- p12 + scale_shape_manual(labels = col_df$group, values = shapes)
  }
  p12 <- p12 + scale_x_continuous(breaks = pretty_breaks(5), expand = expand) 
  p12 <- p12 + scale_y_continuous(breaks = pretty_breaks(5), expand = expand)
  p12 <- p12 + xlab("Tsne Dim1")
  p12 <- p12 + ylab("Tsne Dim2")
  p12 <- p12 + theme_classic() + the
  if(length(unique(egv_order$group)) <= 1){
    p12 <- p12 + theme(legend.position = "none")
  }
  if(border == T){
    p12 <- p12 + theme(panel.border = element_rect(colour = "black", fill=NA, size=border_size))
  }
  pdf(file = paste0(key, ".PC1_PC2.pdf"), width = width, height = height)
  print(p12)
  dev.off()
  png(filename = paste0(key, ".PC1_PC2.png"), width = width, height = height, units = "in", res = 300)
  print(p12)
  dev.off()
  

}
