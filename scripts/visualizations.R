#Circos / diff/meth visual
ideoDMC <- function(obj, chr.len, difference = 25, qvalue = 0.01, title = "test", hyper.col = "magenta", hypo.col = "green") 
{
    myIdeo <- GRanges(seqnames = names(chr.len), ranges = IRanges(start = 1,  width = chr.len))
    seqlevels(myIdeo) = names(chr.len)
    seqlengths(myIdeo) = (chr.len)
    hypo = getMethylDiff(obj, difference = difference, qvalue = qvalue, type = "hypo")
    hyper = getMethylDiff(obj, difference = difference, qvalue = qvalue, type = "hyper")
    g.per = as(hyper, "GRanges")
    seqlevels(g.per, force=TRUE) = seqlevels(myIdeo)
    seqlengths(g.per)=(chr.len)
    values(g.per)$id = "hyper"
    if(nrow(hypo)>0){
      g.po = as(hypo, "GRanges")
      seqlevels(g.po, force=TRUE) = seqlevels(myIdeo)
      seqlengths(g.po)=(chr.len)
      values(g.po)$id = "hypo"
      p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 39, trackWidth = 2)
      p <- p + layout_circle(c(g.po, g.per), geom = "point", size = 1, aes(
        x = midpoint, y = meth.diff, color = id), radius = 25, trackWidth = 30) + scale_colour_manual(values = c(hyper.col, hypo.col))
    }else{
      p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", radius = 29, trackWidth = 2)
      p <- p + layout_circle(g.per, geom = "point", size = 1, aes(
        x = midpoint, y = meth.diff, color = id), radius = 35, trackWidth = 21) + scale_colour_manual(values = hyper.col)
    }
  
    p + layout_circle(myIdeo, geom = "text", aes(
      label = seqnames), vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
}

heatmapMatrixConstruction <- function(mat_meth, plus.minus.all, lookup.tile.ranges, cluster){
 

  #type
  type <- as.character(plus.minus.all)
  type[grepl('1', type)] <- 'Experiment'
  type[grepl('0', type)] <- 'Control'

  #gene_symbol
  gene_symbol <- (as.data.frame(lookup.tile.ranges[2:nrow(lookup.tile.ranges), 2]))

  #generate directions for methylation
  control <- mat_meth[,as.logical(plus.minus.all), drop=F]
  experiment <- mat_meth[,!as.logical(plus.minus.all), drop=F]
  direction = rowMeans(control) - rowMeans(experiment)
  direction = ifelse(direction > 0, "hyper", "hypo")

  #TO BE DETERMINED
  # generate expression matrix
  mat_expr <- NULL
  # matrix for correlation between methylation and expression
  cor_pvalue <- NULL
  # matrix for types of genes
  gene_type <- NULL
  # annotation to genes
  anno_gene <- NULL
  # distance to genes
  dist <- NULL
  # annotation to enhancers
  rand_enhancer <- NULL

  res_list = list()
  res_list$type = type
  res_list$mat_meth = mat_meth
  # res_list$mat_expr = mat_expr
  res_list$direction = direction
  # res_list$cor_pvalue = cor_pvalue
  # res_list$gene_type = gene_type
  # res_list$anno_gene = anno_gene
  # res_list$dist = dist
  # res_list$anno_enhancer = anno_enhancer
  if(cluster==TRUE){
    res_list$column_tree = hclust(dist(t(mat_meth)))
  }
  return (res_list)
}

methylationHeatmap <- function(res_list, title, name, cluster, range){

  ht_global_opt(RESET = TRUE)

  ha = HeatmapAnnotation(df = data.frame(type = res_list$type), 
                         
    col = list(type = c("Experiment" = "pink", "Control" = "royalblue"))
    )

  #MAIN 

  if(cluster==TRUE){
     Heatmap(
        res_list$mat_meth, 
        name = name,
        col = colorRamp2(c(range[1],range[2],range[3]), c("blue", "white", "red")),
        cluster_columns = res_list$column_tree, 
        column_dend_reorder = FALSE, 
        top_annotation = ha, 
        km = 5, 
        column_title = title, 
        column_title_gp = gpar(fontsize = 10)
  ) 
  }else{
      Heatmap(
        res_list$mat_meth, 
        name = name,
        col = colorRamp2(c(range[1],range[2],range[3]), c("blue", "white", "red")),
        km = 5, 
        column_title = title, 
        column_title_gp = gpar(fontsize = 10)
    )
  }

}