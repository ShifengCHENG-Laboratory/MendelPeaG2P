args <- commandArgs(T)
setwd("/vol3/agis/chengshifeng_group/fengcong/zz.my_jupyter_notebook/zhengli_G2B_allele_QTL_NIL/FT-A1-allele-QTL-MTA-NIL")
tree_string <- readLines("TraesCS7A02G115400.g2k.SNPINDEL.AC2.newick",n=1)
tree_string <- readLines(args[1], n = 1)

tree_string <- sub(";+$", "", tree_string)
tree_string <- substr(tree_string, 2, nchar(tree_seq)-1)
# cat(tree_string)

library(stringr)

substr_count <- function(string, pattern) {
  return(sum(str_count(string, fixed(pattern))))
}


find_comma_index <- function(tree_seq) {
  comma_index_list <- c()
  for (index in 1:nchar(tree_seq)) {
    if (substr(tree_seq, index, index) == ",") {
      comma_index_list <- c(comma_index_list, index)
    }
  }
  
  i_list <- c()
  for (i in comma_index_list) {

    if (substr_count(substr(tree_seq, 1, i-1), "(") == substr_count(substr(tree_seq, 1, i-1), ")") &&
        substr_count(substr(tree_seq, i + 1, nchar(tree_seq)), ")") == substr_count(substr(tree_seq, i + 1, nchar(tree_seq)), "(")) {
      i_list <- c(i_list, i)
    }
  }
  
  return(i_list)
}



rm_bracket2 <- function(tree_seq) {
  if (substr(tree_seq, 1, 1) == "(") {
    pos <- str_locate_all(tree_seq, fixed(")"))[[1]]
    # 取出最后一个右括号的位置
    rbound <- tail(pos, 1)[1]
    if (rbound == nchar(tree_seq)) {
      return(list(substr(tree_seq, 2, nchar(tree_seq) - 1), 1))
    } else {
      return(list(substr(tree_seq, 2, rbound-1), as.numeric( strsplit(tree_seq, ":")[[1]][length(strsplit(tree_seq, ":")[[1]])] )))
    }
  } else {
    if (grepl(":", tree_seq)) {
      return(list(tree_seq, as.numeric(strsplit(tree_seq, ":")[[1]][length(strsplit(tree_seq, ":")[[1]])])))
    } else {
      return(list(tree_seq, 1))
    }
  }
}

post_traverse_tree <- function(tree_seq, cur_branch_len) {
  comma_index <- find_comma_index(tree_seq)
  comma_index_len <- length(comma_index)

  if (comma_index_len == 0) {
    sample_name <- strsplit(tree_seq, ":")[[1]][1]
    return(list(list(sample_name), cur_branch_len))
  }

  seq_list <- list()
  pre_start <- 1
  for (i in 1:comma_index_len) {
    indx <- comma_index[i]
    sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, indx - 1))
    seq_list[[i]] <- sub_tree_seq_branch_len
    pre_start <- indx + 1
    if (i == comma_index_len) {
      sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, nchar(tree_seq)))
      seq_list[[i + 1]] <- sub_tree_seq_branch_len
    }
  }

  sample_list <- list()
  cur_max_branch_len <- 0
  for (i in 1:length(seq_list)) {
    sub_tree_seq <- seq_list[[i]][[1]]
    sub_branch_len <- seq_list[[i]][[2]]
    sub_sample_list_max_depth <- post_traverse_tree(sub_tree_seq, cur_branch_len + sub_branch_len)
    sub_sample_list <- sub_sample_list_max_depth[[1]]
    sub_max_depth <- sub_sample_list_max_depth[[2]]
    sample_list <- c(sample_list, sub_sample_list)
    cur_max_branch_len <- max(cur_max_branch_len, sub_max_depth)
  }

  return(list(sample_list, cur_max_branch_len))
}

library(graphics)
draw_line <- function(x1, y1, x2, y2, lty = 1) {
  if (x1 == x2 && y1 == y2) {
    return()
  }
  print(paste(x1, y1, x2, y2, sep = ","))
  lines(x = c(x1, x2), y = c(y1, y2), lty = lty)
}



post_traverse_tree_and_plot <- function(tree_seq, samples_up_to_you, father_depth, cur_depth, cell_unit) {
  comma_index <- find_comma_index(tree_seq)
  comma_index_len <- length(comma_index)
  
  if (comma_index_len == 0) {
    sample_name <- strsplit(tree_seq, ":")[[1]][1]
    samples_up_to_you <- append(samples_up_to_you, sample_name)
    row_num <- length(samples_up_to_you)
    father_col <- as.integer(father_depth / cell_unit)
    cur_col <- as.integer(cur_depth / cell_unit)
    
    draw_line(father_col , -row_num, cur_col , -row_num)
    # print(paste0("father_col: ", father_col, " cur_col: ", cur_col, " row_num: ", row_num))
    # print(samples_up_to_you)
    return(list(samples_up_to_you, row_num))
  }
  
  seq_list <- list()
  pre_start <- 1
  for (i in 1:comma_index_len) {
    indx <- comma_index[i]
    sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, indx - 1))
    seq_list[[i]] <- sub_tree_seq_branch_len
    pre_start <- indx + 1
    if (i == comma_index_len) {
      sub_tree_seq_branch_len <- rm_bracket2(substr(tree_seq, pre_start, nchar(tree_seq)))
      seq_list[[i + 1]] <- sub_tree_seq_branch_len
    }
  }
  
  
  row_num_list <- list()
  for (i in 1:length(seq_list)) {
    sub_tree_seq <- seq_list[[i]][[1]]
    sub_branch_len <- seq_list[[i]][[2]]
    sub_sample_list_row_num <- post_traverse_tree_and_plot(sub_tree_seq, samples_up_to_you, cur_depth, cur_depth + sub_branch_len, cell_unit)
    samples_up_to_you <- sub_sample_list_row_num[[1]]
    row_num <- sub_sample_list_row_num[[2]]
    row_num_list[[i]] <- row_num
  }
  
  father_col <- as.integer(father_depth / cell_unit)
  cur_col <- as.integer(cur_depth / cell_unit)
#   col_num_list <- (father_col + 1):(cur_col + 1)
  row_num <- (row_num_list[[1]] + row_num_list[[length(row_num_list)]]) / 2

  # plot node
  #    |
  #----|
  #    |

  draw_line(cur_col, -row_num_list[[1]], cur_col, -row_num_list[[length(row_num_list)]])
  draw_line(father_col, -row_num, cur_col, -row_num)
  

  #plot root
  # |
  # |
  # |
#   if (father_depth == 0 && cur_depth == 0) {
#     draw_line(1, -row_num_list[[1]], 1, -row_num_list[[length(row_num_list)]])
#   }
  return(list(samples_up_to_you, row_num))
}




tree_seq = tree_string
# Example usage
tree_seq <- "(K:0.1,B:0.2):0.3,(C:0.3,D:0.4):0.5"
sample_list_max_depth <- post_traverse_tree(tree_seq, 0)
sample_list <- sample_list_max_depth[[1]]
max_depth <- sample_list_max_depth[[2]]

print(sample_list)
print(max_depth)

tree_width=1000
cell_unit = (max_depth*1.0)/(tree_width)

pdf("tree.pdf",height = 11.7,width = 8.3)
par(plt = c(0, 1, 0, 1))
plot.new()
plot.window(xlim = c(0, 1001), ylim = c(-1047, 0))
# plot(1,1,xlim=c(0,1001),ylim=c(1 ,-6),type="n",axes=F,xlab="",ylab="")
# 绘制树结构
result <- post_traverse_tree_and_plot(tree_seq, list(), 0, 0, cell_unit)
dev.off()
# 显示绘制的图形
draw_line(0,0,1001,0)
draw_line(0,0,0,-6)

box()
title("Phylogenetic Tree")
# substr_count(tree_string, ")")


