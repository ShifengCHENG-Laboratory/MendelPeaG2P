gene_structure <- function(geneStructure_df){
    # geneStructure file(CDS,UTR,intron,upstream,downstream):
    # start     end     feature
    # 1         50      CDS
    # ...

    # markPos file: 
    # position
    # 49

    ## get max and min position for this gene
    min_pos = min(geneStructure_df$start)
    max_pos = max(geneStructure_df$end)
    region_length = max_pos - min_pos+1

    ##setting height(y)
    height_UTR = 0.2
    height_intron = 0.07
    height_CDS = 0.5
    height_other = 0.01

    base_y = 20

    for(i in seq(1,nrow(geneStructure_df),1)){
        plotx = c()
        ploty = c()
        cl = NA
        if (geneStructure_df$feature[i] == "CDS"){
            ploty = c(base_y+height_CDS, base_y+height_CDS , base_y-height_CDS ,base_y-height_CDS )
        }else if(geneStructure_df$feature[i] == "UTR"){
            ploty = c(base_y+height_UTR, base_y+height_UTR , base_y-height_UTR ,base_y-height_UTR )
        }else if(geneStructure_df$feature[i] == "intron"){
            ploty = c(base_y+height_intron, base_y+height_intron , base_y-height_intron ,base_y-height_intron )
            cl = "black"
        }else{
            ploty = c(base_y+height_other, base_y+height_other , base_y-height_other ,base_y-height_other )
        }
        
        plotx = c(geneStructure_df$start[i]-min_pos, geneStructure_df$end[i]-min_pos, geneStructure_df$end[i]-min_pos, geneStructure_df$start[i]-min_pos)
        polygon(plotx,ploty,col = cl,lwd=1)
    }



    ## plot axis
    lines(c(0,region_length), c(base_y-1, base_y-1),col = "black" ,lwd=1)
    for(i in seq(0,region_length,500)){
      lines(c(i,i),c(base_y-1,base_y-1.2),lwd=1)
      text(i,base_y-1.5,i+min_pos,cex=0.5)
    }

}

args = commandArgs(T)
#genestructure , markpos , ouput
gene_structure(args[1],args[2],args[3])