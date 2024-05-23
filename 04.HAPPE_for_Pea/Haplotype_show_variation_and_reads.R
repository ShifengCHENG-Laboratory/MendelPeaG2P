#!/usr/bin/env Rscript
library(optparse)
library(graphics)
library(stringr)
paper_size_height =  6
paper_size_width = 4 #11.7
phe_colors <- c("#FFC000","#029676","#8B4513","#AD84C6","#FB8072", "#D092A7","#0989B1","#809EC2", "#A5B592")

main_function <- function(var_file, sample_file, gff_file,  region, phe_file, is_quantity, output_prefix) {
    # 在这里编写您的代码逻辑
    # cat("heat_file:", heat_file, "\n")
    # cat("tree_file:", tree_file, "\n")
    # cat("cluster_file:", cluster_file, "\n")
    # cat("gff_file:", gff_file, "\n")
    # cat("sample_color_file:", sample_color_file, "\n")
    # cat("region:", region, "\n")
    # cat("output_prefix:", output_prefix, "\n")

    ##### 处理var_file #####
    ## 定义var file里面的type 与 color的对应关系
    ## S：#FF0000，I：#FFFF00，D：#0000FF，Q：#808080
    var_color_mapping <- list(
        "S" = "#FF0000",
        "D" = "#ADD8E6" , #"I" = "#FFFF00", 
        "I" = "#0000FF",
        "Q" = "#CCCCCC",
        "E" = "#FFA500"
    )

    var_data = process_var_file(var_file)


    ##### 处理sample file #####
    sample_data <- read.table(sample_file, header = F, stringsAsFactors = FALSE,check.names = FALSE) ## 他只有一列，转为vector吧
    colnames(sample_data) <- c("sample")
    sample_list <- as.vector(sample_data$sample)


    ##### 处理 region #####
    # interval_size = get_interval_size(region)
    # 分割字符串并提取区间部分
    interval_string <- unlist(strsplit(region, ":"))[2]
    # 分割区间字符串并转换为数值
    interval_bounds <- as.numeric(unlist(strsplit(interval_string, "-")))
    # 计算区间大小（包含边界）
    interval_size <- interval_bounds[2] - interval_bounds[1] + 1
    interval_start = interval_bounds[1]
    interval_end = interval_bounds[2]
    interval_vector = c(interval_start, interval_end, interval_size)

    ## 计算各个部分占的长度
    left_white_space = interval_size * 0
    var_phenotype_space = interval_size * 0.01
    phenotype_space = interval_size * 0.01
    right_white_space = interval_size * 0

    max_x = interval_size + left_white_space + var_phenotype_space + phenotype_space + right_white_space
    white_space <<- nrow(sample_data) *(1/21)
    max_y = nrow(sample_data) + white_space
    
    #####  处理gff file #####
    combined_data = process_gff_file3(gff_file, region)


    ##### 准备画图 #####
    pdf(file=paste(output_prefix,"pdf",sep="."), height = paper_size_height, width = paper_size_width)
    par(mar=c(1, 1, 1, 1))
    # plot(1,1,xlim=c(0,max_x),ylim=c(-max_y, 0),type="n",axes=F,xlab="",ylab="")
    plot(1, 1, xlim=c(-left_white_space, max_x-left_white_space),ylim=c(-nrow(sample_data), white_space), type = "n", axes = FALSE, xlab = "", ylab = "")
    # -----------> x = c(0,max_x)
    # |
    # |
    # ↓
    # y = c(0,-max_y)

    ###### 画 每个样本
    draw_var_matrix(var_data, interval_vector, var_color_mapping, sample_data)

    ### 画基因结构
    draw_gene_structure3(combined_data)

    # 检查是否提供了 --phenotype 和 --is_quantity 选项
    #### 画表型数据 ####
    if (!is.null(phe_file) && !is.null(is_quantity)) {
      phe_data <- read.table(phe_file, header = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
      names(phe_data) = c("sample", "phe")
      max_phe <- max(phe_data$phe, na.rm = TRUE)
      phe_data_dict <- setNames(phe_data$phe, phe_data$sample)
      draw_trait_phenotype(phe_data_dict, interval_size, sample_list, is_quantity, max_phe,var_phenotype_space, phenotype_space)
      
    }

    dev.off()
}

draw_trait_phenotype <- function(phe_data_dict, interval_size, sample_list, is_quantity, max_phe,var_phenotype_space,phenotype_space) {
  n_rows <- length(sample_list)
  rect_height <- 0.8
  # 遍历样本绘制矩形
  if (is_quantity) {
    phe_color <- phe_colors[4]
    for (i in 1:n_rows) {
      # 获取单元格的颜色
      sample_name <- sample_list[i]
      # sample_phe <- phe_data_dict[[sample_name]]
      if (sample_name %in% names(phe_data_dict)) {
        sample_phe <- phe_data_dict[[sample_name]]
      } else {
        sample_phe <- NA
      }
      if (!is.na(sample_phe)) {
        # 计算矩形的边界
        x_left <- interval_size + var_phenotype_space
        x_right <- x_left + ((phenotype_space * sample_phe) / max_phe)
        y_top <- -1 * i + rect_height/2
        y_bottom <- y_top - rect_height
        
        # 绘制无边框矩形
        rect(x_left, y_bottom, x_right, y_top, col = phe_color, border = NA)
        # text(x_right + text_width * 1.15, (y_top + y_bottom) / 2, round(sample_phe, 2), cex = 0.16, pos = 4)
      }
    }
  } else {
    for (i in 1:n_rows) {
      sample_name <- sample_list[i]
      if (sample_name %in% names(phe_data_dict)) {
        sample_phe <- phe_data_dict[[sample_name]]
      } else {
        sample_phe <- NA
      }
      if (!is.na(sample_phe)) {
        # 获取单元格的颜色
        phe_color <- ifelse(as.numeric(sample_phe) == 0, "#000000", phe_colors[as.numeric(sample_phe)])
        # 计算矩形的边界
        x_left <- interval_size + var_phenotype_space
        x_right <- x_left + phenotype_space / 2
        y_top <- -1 * i + rect_height/2 
        y_bottom <- y_top - rect_height
        
        # 绘制无边框矩形
        rect(x_left, y_bottom, x_right, y_top, col = phe_color, border = NA)
      }
    }
  }
}

draw_gene_structure3 <- function(combined_data) {
  # Prepare data
  gene_structure_df <- combined_data
  names(gene_structure_df)[4] <- "feature"
  #   print(combined_data)
  gene_height = white_space
  # Set parameters
  height_UTR <- 5#gene_height / 10
  height_intron <- 1#gene_height / 25
  height_CDS <- 5#gene_height / 10
  #   height_other <- 1
  base_y <- gene_height / 2
  arrow_y <- base_y
  
  min_pos <- min(gene_structure_df$start) 
  base_x = 0
  # Draw gene structure
  for (i in seq(1, nrow(gene_structure_df), 1)) {
    arrow_y <- base_y
    plotx <- c()
    ploty <- c()
    cl <- NA
    if (gene_structure_df$feature[i] == "CDS") {
      ploty <- c(base_y + height_CDS, base_y + height_CDS, base_y - height_CDS, base_y - height_CDS)
      cl <- "#007BA7"
      arrow_y <- base_y - height_CDS/2
    } else if (gene_structure_df$feature[i] == "UTR") {
      ploty <- c(base_y + height_UTR, base_y + height_UTR, base_y - height_UTR, base_y - height_UTR)
      arrow_y <- base_y - height_UTR/2
    } else if (gene_structure_df$feature[i] == "intron") {
      ploty <- c(base_y + height_intron, base_y + height_intron, base_y - height_intron, base_y - height_intron)
      cl <- "black"
      arrow_y <- base_y - height_intron/2
    } else if (gene_structure_df$feature[i] == "intergenic") {
      ploty <- c(base_y, base_y)
    } else { ## gene name
      gene_name <- gene_structure_df$feature[i]
      gene_start <- gene_structure_df$start[gene_structure_df$feature == gene_name]
      gene_end <- gene_structure_df$end[gene_structure_df$feature == gene_name]
      gene_center <- (gene_start + gene_end) / 2 - min_pos
      gene_name = paste(gene_name,"(",gene_structure_df$strand[i],")",sep="")
      text(base_x+gene_center, base_y + gene_height*(5/5), gene_name, cex = 0.7)
      next
    }

    
    plotx <- c(base_x+gene_structure_df$start[i] - min_pos, base_x+gene_structure_df$end[i] - min_pos, base_x+gene_structure_df$end[i] - min_pos, base_x+gene_structure_df$start[i] - min_pos)
    
    if (gene_structure_df$feature[i] == "intergenic") {
      #  lines(plotx[1:2], ploty[1:2], col = "black", lwd = 1) 
      # 暂时不画intergenic
      
    } else if (gene_structure_df$feature[i] == "intron"){
        ## 如果是intron 画一个曲线
        intro_len = gene_structure_df$end[i] - gene_structure_df$start[i]
        lines(c(base_x+gene_structure_df$start[i] - min_pos, base_x+gene_structure_df$start[i] - min_pos + intro_len/2), c(base_y,base_y-height_intron*(1)), col = "black", lwd = 1)
        lines(c(base_x+gene_structure_df$start[i] - min_pos + intro_len/2, base_x+gene_structure_df$end[i] - min_pos), c(base_y-height_intron*(1),base_y), col = "black", lwd = 1)

    } else {
      polygon(plotx, ploty, col = cl, lwd = 0.5)
      arrow_x <- (gene_structure_df$start[i] + gene_structure_df$end[i]) / 2 - min_pos
      if (gene_structure_df$strand[i] == "+") {
        arrow_direction <- ">"
      } else {
        arrow_direction <- "<"
      }
      # 用红色标识方向
    #   text(base_x + arrow_x, arrow_y, arrow_direction, cex = 0.7, col = "red")
    }
  }

#   # Draw axis
#   region_length <- max(gene_structure_df$end) - min_pos + 1
#   lines(c(base_x+0, base_x+region_length), c(base_y + 10, base_y + 10), col = "black", lwd = 1)
  
#   ## 确定一个合适的间隔，stepsize，这个值是1000的倍数，保证这个区间最多只有10个刻度
#   stepsize = 1000
#   while (region_length/stepsize > 10){
#       stepsize = stepsize + 1000
#   }


#   for (i in seq(0, region_length, stepsize)) {
#     lines(c(base_x+i, base_x+i), c(base_y + 10, base_y + 11), lwd = 1)
#     x = formatC(i + min_pos, format = "f", big.mark = ",", digits = 0)
#     # print(x)
#     text(base_x+i, base_y + 15, x, cex = 0.5)
#   }

}

draw_var_matrix <- function(var_data, interval_vector, var_color_mapping, sample_data){
    interval_start = interval_vector[1]
    interval_end = interval_vector[2]
    ## 遍历每一个样本
    for (i in 1:nrow(sample_data)){
        sample = sample_data[i, "sample"]
        sample_var_data = var_data[var_data$sample == sample,]
        # rect(1, -i-0.05, interval_vector[3], -i+0.05, col="#CCCCCC", border=NA)

        # print(sample)
        # print(sample_var_data)
        ## 画一个样本的变异，遍历sample_var_data,要优先画var 为Q的变异
        for (j in 1:nrow(sample_var_data)){
            
            var = as.character(sample_var_data[j, "type"])
            if (var != "Q"){
                next
            }
            start = as.numeric(sample_var_data[j, "start"] - interval_start)
            end = as.numeric(sample_var_data[j, "end"] - interval_start)
            color = as.character(var_color_mapping[var])

            rect(start, -i-0.4, end, -i+0.4, col=color, border=NA)

        }
        for (j in 1:nrow(sample_var_data)){
            
            var = as.character(sample_var_data[j, "type"])
            if (var == "Q"){
                next
            }
            start = as.numeric(sample_var_data[j, "start"] - interval_start)
            end = as.numeric(sample_var_data[j, "end"] - interval_start)
            color = as.character(var_color_mapping[var])

            rect(start-1, -i-0.4, end+1, -i+0.4, col=color, border=NA)
            ## 如果var 为I 需要在上下各画一个横着的矩形
            if (var == "I"){
                rect(start-2, -i-0.4, end+2, -i-0.39, col=color, border=NA)
                rect(start-2, -i+0.39, end+2, -i+0.4, col=color, border=NA)
            }
        }

    }
}




process_var_file <- function(var_file) {
  var_data <- read.table(var_file, header = F, stringsAsFactors = FALSE,check.names = FALSE)
  colnames(var_data) <- c("sample", "type", "start", "end", "homohet")
  return(var_data)
}

get_interval_size <- function(region){
    # 分割字符串并提取区间部分
    interval_string <- unlist(strsplit(region, ":"))[2]
    # 分割区间字符串并转换为数值
    interval_bounds <- as.numeric(unlist(strsplit(interval_string, "-")))
    # 计算区间大小（包含边界）
    interval_size <- interval_bounds[2] - interval_bounds[1] + 1
    # cat("interval_size:", interval_size, "\n")

    return(interval_size)
}

### 仅仅适用于豌豆 gff
process_gff_file3 <- function(gff_file, region){
  gff_data <- read.table(gff_file, sep="\t", header = FALSE, stringsAsFactors = FALSE,
    quote = "", comment.char = "#",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))

  # 筛选 region 对应的染色体数据
  # 潜在bug,提供的区间 和 gff 文件中的染色体名称不一致
  region_seqid <- unlist(strsplit(region, ":"))[1]
  gff_data_filtered <- gff_data[gff_data$seqid == region_seqid,]

  # 筛选 region 对应的区间数据
  interval_string <- unlist(strsplit(region, ":"))[2]
  interval_bounds <- as.numeric(unlist(strsplit(interval_string, "-")))
  gff_data_filtered <- gff_data_filtered[gff_data_filtered$start >= interval_bounds[1] & gff_data_filtered$end <= interval_bounds[2],]

  if (nrow(gff_data_filtered) > 0) {
    # 提取 gene、CDS、UTR 和 exon 数据
    gene_data <- gff_data_filtered[gff_data_filtered$type == "gene", c("start", "end", "attributes", "strand")]
    cds_data <- gff_data_filtered[gff_data_filtered$type == "CDS", c("start", "end", "strand", "attributes")]
    utr_data <- gff_data_filtered[gff_data_filtered$type %in% c("five_prime_UTR", "three_prime_UTR"), c("start", "end", "strand", "attributes")]
    exon_data <- gff_data_filtered[gff_data_filtered$type == "exon", c("start", "end", "strand", "attributes")]
    mrna_data <- gff_data_filtered[gff_data_filtered$type == "mRNA", c("start", "end", "strand", "attributes")]
    # 获取第一个转录本的ID
    # 潜在bug: 并不是所有物种都是用这种正则能匹配到 gene id
    first_transcript_ids <- sapply(strsplit(mrna_data$attributes, ";"), function(x) {
      id_field <- x[grep("^ID=", x)]
      if (length(id_field) > 0) {
        sub("^ID=rna-gnl\\|WGS:JAMSHJ\\|", "", id_field)
      } else {
        NA
      }
    })[1]

    print(first_transcript_ids)
    # gene_data <- subset(gene_data, select = -attributes)
    cds_data_filtered <- data.frame()
    utr_data_filtered <- data.frame()
    exon_data_filtered <- data.frame()
    intron_data <- data.frame()
    for (transcript_id in first_transcript_ids) {
      cds_data_filtered <- rbind(cds_data_filtered,
                                gff_data_filtered[gff_data_filtered$type == "CDS" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      utr_data_filtered <- rbind(utr_data_filtered,
                                gff_data_filtered[gff_data_filtered$type %in% c("five_prime_UTR", "three_prime_UTR") & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      exon_data_filtered <- rbind(exon_data_filtered,
                                gff_data_filtered[gff_data_filtered$type == "exon" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")])
      # 计算 intron 数据
      tmp_exon_data = gff_data_filtered[gff_data_filtered$type == "exon" & grepl(transcript_id, gff_data_filtered$attributes), c("start", "end", "strand")]
      if (nrow(tmp_exon_data) > 1) {
        for (i in 1:(nrow(tmp_exon_data) - 1)) {
          if (tmp_exon_data[i, "end"] + 1 < tmp_exon_data[i + 1, "start"]) {
            intron_data <- rbind(intron_data, data.frame(start = tmp_exon_data[i, "end"] + 1, end = tmp_exon_data[i + 1, "start"] - 1))
          } else {
            intron_data <- rbind(intron_data, data.frame(start = tmp_exon_data[i, "start"] + 1, end = tmp_exon_data[i + 1, "end"] - 1))
          }
        }
      }
    }
    # print(intron_data)
    # 用筛选后的数据更新CDS、UTR和exon数据
    cds_data <- cds_data_filtered
    utr_data <- utr_data_filtered
    exon_data <- exon_data_filtered
    # # 筛选第一个转录本的CDS、UTR和exon数据
    # if (nrow(cds_data) > 0) {
    #   cds_data <- cds_data[grep(paste0("Parent=", first_transcript_id), cds_data$attributes), ]
    # }
    # if (nrow(utr_data) > 0) {
    #   utr_data <- utr_data[grep(paste0("Parent=", first_transcript_id), utr_data$attributes), ]
    # }
    # if (nrow(exon_data) > 0) {
    #   exon_data <- exon_data[grep(paste0("Parent=", first_transcript_id), exon_data$attributes), ]
    # }

    # 删除不再需要的attributes列
    # cds_data <- subset(cds_data, select = -attributes)
    # utr_data <- subset(utr_data, select = -attributes)
    # exon_data <- subset(exon_data, select = -attributes)

    # 从 attributes 列中提取 ID 并分配给 gene_data 的 type 列
    if (nrow(gene_data) > 0) {
        gene_data$type <- sapply(strsplit(gene_data$attributes, ";"), function(x) {
        id_field <- x[grep("^ID=", x)]
        ## gene id 加上.1 
        id_field <- paste0(id_field, ".1")
        if (length(id_field) > 0) {
          sub("^ID=", "", id_field)
        } else {
          NA
        }
      })
      gene_data <- subset(gene_data, select = -attributes)
    }

    

    # 为所有数据添加 type 列，前提是数据框不为空
    if (nrow(cds_data) > 0) cds_data$type <- "CDS"
    if (nrow(utr_data) > 0) utr_data$type <- "UTR"
    if (nrow(exon_data) > 0) exon_data$type <- "exon"
    if (nrow(intron_data) > 0) intron_data$type <- "intron"

    ## 给intron_data添加strand列，根据其他data的strand列来判断
    if (nrow(intron_data) > 0) {
      intron_data$strand <- exon_data$strand[1]
      # 并且将列排一下序，让其和其他data的列顺序一致
      intron_data <- intron_data[,c("start", "end", "strand", "type")]
    }

    # 合并数据
    combined_data <- rbind(gene_data, cds_data, utr_data, intron_data)

    # 添加 intergenic 区域
    # 潜在bug： gene区间和当前转录本占用的区间不一定一致，因为gene区间是所有转录本的区间的并集
    # 优化参考：https://en.define.sh/How-to-Calculate-Intersection-Union-And-Difference-of-Intervals
    intergenic_start <- interval_bounds[1]
    intergenic_end <- interval_bounds[2]
    intergenic_data <- data.frame()
    if (nrow(gene_data) > 0) {
      for (i in 1:nrow(gene_data)) {
        if (gene_data[i, "start"] > intergenic_start) {
          intergenic_data <- rbind(intergenic_data, data.frame(start = intergenic_start, end = gene_data[i, "start"] - 1, type = "intergenic"))
        }
        intergenic_start <- gene_data[i, "end"] + 1
      }
    }
    if (intergenic_start <= intergenic_end) {
      intergenic_data <- rbind(intergenic_data, data.frame(start = intergenic_start, end = intergenic_end, type = "intergenic"))
    }

    # 合并 intergenic_data 并按照 start 列排序
    if (nrow(intergenic_data) > 0) {
        intergenic_data$strand <- "x"
        combined_data <- rbind(combined_data, intergenic_data)
    }
    combined_data <- combined_data[order(combined_data$start), ]
  } else {
    combined_data <- data.frame(start = interval_bounds[1], end = interval_bounds[2], strand <- "x", type = "intergenic")
  }
  return(combined_data)
}



# 处理命令行参数
process_command_line <- function() {
    # 定义命令行选项
    option_list <- list(
        make_option(c("-H", "--var_file"), type = "character", help = "Path to the var_file."),
        make_option(c("-s", "--sample_file"), type = "character", help = "Path to the sample_file."),
        make_option(c("-g", "--gff_file"), type = "character", help = "Path to the gff_file."),
        make_option(c("-r", "--region"), type = "character", help = "Region."),
        make_option(c("-p", "--phenotype"), type = "character", default = NULL, help = "Phenotype data."),
        make_option(c("-f", "--is_quantity"), type = "logical", default = NULL, help = "Is it a quantitative trait?"),
        make_option(c("-o", "--output_prefix"), type = "character", help = "Output file prefix.")
    )

    # 解析命令行参数
    args <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)

    if (!is.null(args$help) && args$help) {
        print_help(OptionParser(option_list = option_list))
        quit(save = "no", status = 0)
    }

    # 检查参数是否完整
    required_args <- c("var_file", "sample_file", "gff_file", "region", "output_prefix")
    for (arg_name in required_args) {
        if (is.null(args$options[[arg_name]])) {
            cat("Error: Argument", arg_name, "is missing.\n")
            print_help(OptionParser(option_list = option_list))
            quit(save = "no", status = 1)
        }
        if ( !(arg_name %in% c("region", "output_prefix")) && !file.exists(args$options[[arg_name]])) {
            cat("Error: File", args$options[[arg_name]], "for argument", arg_name, "does not exist.\n")
            print_help(OptionParser(option_list = option_list))
            quit(save = "no", status = 1)
        }
    }
    if ((!is.null(args$options$phenotype) && is.null(args$options$is_quantity)) || is.null(args$options$phenotype) && !is.null(args$options$is_quantity)) {
        cat("The phenotype file and the is_quantity value must both exist.\n")
        print_help(OptionParser(option_list = option_list))
        quit(save = "no", status = 1)
    }
    main_function(args$options$var_file, args$options$sample_file,  args$options$gff_file,
                  args$options$region, args$options$phenotype,
                  args$options$is_quantity, args$options$output_prefix)
}

# 如果脚本是通过 Rscript 直接运行的，处理命令行参数并调用主函数
if (interactive() == FALSE) {
    process_command_line()
}