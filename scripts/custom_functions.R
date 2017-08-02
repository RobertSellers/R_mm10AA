annotateRefSymbol <- function(x, gene.parts, cpg.obj, pct){

    #windows
    gRange <- as(x,"GRanges")
    gRange.annot = genomation::annotateWithGeneParts(gRange, gene.parts, intersect.chr = TRUE)
    genomation::getTargetAnnotationStats(gRange.annot, percentage = TRUE, precedence = TRUE)
    genomation::plotTargetAnnotation(gRange.annot, precedence = TRUE, main = paste("Differential Methylation Annotation", pct, "%"))

    #cpg islands
    gRange.cpg <- genomation::annotateWithFeatureFlank(gRange, cpg.obj$CpGi, cpg.obj$shores, feature.name = "CpGi", flank.name = "shores")
    genomation::getTargetAnnotationStats(gRange.cpg, percentage = TRUE, precedence = TRUE)
    genomation::plotTargetAnnotation(gRange.cpg, precedence = TRUE, main = paste("Differential Methylation CpGi", pct, "%"))

}

filterAndNormalize <- function(x){
    filename <- gsub("[.]","",deparse(substitute(x))) 
    destfile <- paste("./save/",filename,"_filt_norm.rds",sep="")
    if(!file.exists(destfile)){
        filt <- filterByCoverage(x, lo.count=10, lo.perc=NULL, hi.count=NULL, hi.perc=99.9)
        filt.norm = normalizeCoverage(filt, method = "median")
        saveRDS(filt.norm, destfile)
    }else{
        filt.norm <- readRDS(destfile)
    }
    return (filt.norm)
}

mergeAndMethylate <- function(x){
    filename <- gsub("[.]","_",deparse(substitute(x))) 
    destfile <- paste("./save/",filename,"_meth.rds",sep="")
    if(!file.exists(destfile)){
        meth <- methylKit::unite(x, destrand=FALSE)
        saveRDS(meth, destfile)
    }else{
        meth<-readRDS(destfile)
    }
    return (meth)
}

cPercentages <- function(all.filt.norm.meth,lookup.tile.ranges){
    destfile <- "cPercentages.rds"
    if(!file.exists(destfile)){

        # x.pct.c <- all.filt.norm.meth %>%
            # group_by(chr, start) %>%
            # group_by(gr=cut(start, breaks= geneRange, right=F)) %>%
            # mutate(pctC_1 = numCs1/coverage1) %>%
            # mutate(pctC_2 = numCs2/coverage2) %>%
            # mutate(pctC_3 = numCs3/coverage3) %>%
            # mutate(pctC_4 = numCs4/coverage4) %>%
            # mutate(pctC_5 = numCs5/coverage5) %>%
            # mutate(pctC_6 = numCs6/coverage6) %>%
            # mutate(pctC_7 = numCs7/coverage7) %>%
            # mutate(pctC_8 = numCs8/coverage8) %>%
            # mutate(pctC_9 = numCs9/coverage9) %>%
            # mutate(pctC_10 = numCs10/coverage10) %>%
            # mutate(pctC_11 = numCs11/coverage11) %>%
            # mutate(pctC_12 = numCs12/coverage12)
            saveRDS(x.pct.c, destfile)
        }else{
            x.pct.c<-readRDS(destfile)
        }
    return (x.pct)
}