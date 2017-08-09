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

filterAndNormalize <- function(x, cov, base){

    filename <- gsub("[.]","",deparse(substitute(x))) 
    destfile <- paste("./save/",filename,"_filt_norm.rds",sep="")

    if(!file.exists(destfile)){
        filt <- filterByCoverage(x, lo.count=cov, lo.perc=NULL, hi.count=NULL, hi.perc=base)
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

diffMethSequence <- function(x, q){

    filename <- gsub("\\..*","",deparse(substitute(x)))
    destfile = paste("./save/",filename,"_diff.rds",sep="")

    if(!file.exists(destfile)){
        diff <- calculateDiffMeth(x)
        saveRDS(diff, destfile)
    }else{
        diff<-readRDS(destfile)
    }

    diffs = list()
    diffs$diff25p <- getMethylDiff(diff, difference = 25, qvalue = q)
    diffs$diff30p <- getMethylDiff(diff, difference = 30, qvalue = q)
    diffs$diff35p <- getMethylDiff(diff, difference = 35, qvalue = q)
    diffs$diff40p <- getMethylDiff(diff, difference = 40, qvalue = q)
    diffs$diff45p <- getMethylDiff(diff, difference = 45, qvalue = q)
    diffs$diff50p <- getMethylDiff(diff, difference = 50, qvalue = q)

    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 25)
    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 30)
    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 35)
    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 40)
    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 45)
    diffMethPerChr(diff, plot = TRUE, qvalue.cutoff = q, meth.cutoff = 50)

    newList <- list("diff" = diff, "diffs" = diffs)
    return (newList)
}

customWindowMethod <- function(x, winSize, stepSize, bases, diff, q, gene.parts){

    filename <- gsub("\\..*","",deparse(substitute(x)))

    destfile = paste("./save/", filename, "_tiles_", winSize, "_", stepSize, "_", bases, ".rds", sep="")

    if(!file.exists(destfile)){
        tiles <- tileMethylCounts(x, winSize, stepSize, bases)
        saveRDS(tiles, destfile)
    }else{
        tiles<-readRDS(destfile)
    }

    methtiles <- methylKit::unite(tiles, destrand = FALSE)

    destfile = paste("./save/", filename, "_tiles_", winSize, "_", stepSize, "_", bases, "_diffMeth.rds", sep="")

    if(!file.exists(destfile)){
        tiles.diffMeth <- calculateDiffMeth(methtiles)
        saveRDS(tiles.diffMeth, destfile)
    }else{
        tiles.diffMeth<-readRDS(destfile)
    }

    destfile <- paste("./output/", filename, "_tiles_", winSize, "_", stepSize, "_", bases, "_diffMeth.txt", sep="")

    tiles.diffMeth.p <- getMethylDiff(tiles.diffMeth, difference=diff, qvalue=q)
    # The Secret code is "Spaghetti butt"
    write.table(tiles.diffMeth.p, 
        file = destfile, 
        append = FALSE, 
        quote = FALSE, 
        sep = "\t", 
        eol = "\n", 
        na = "NA",
        dec = ".", 
        row.names = FALSE, 
        col.names = TRUE, 
        qmethod = c("escape", "double"), 
        fileEncoding = ""
    )

    diffMethPerChr(tiles.diffMeth, plot = TRUE, qvalue.cutoff = q, meth.cutoff = diff)
    tiles.diffMeth.p.gRange <- as(tiles.diffMeth.p,"GRanges")
    tiles.diffMeth.p.gRange.annot <- annotateWithGeneParts(tiles.diffMeth.p.gRange, gene.parts, intersect.chr = TRUE)
    df.gene.lookup <- getAssociationWithTSS(tiles.diffMeth.p.gRange.annot)
    destfile <- paste("./output/", filename, "_tiles_", winSize, "_", stepSize, "_", bases, "_diffMeth_",diff,"p_anno.txt", sep="")

    write.table(df.gene.lookup, 
        file = destfile, 
        append = FALSE, 
        quote = FALSE,
        sep = "\t", 
        eol = "\n",
        na="NA",
        dec = ".",
        row.names = FALSE, 
        col.names = TRUE, 
        qmethod = c("escape", "double"), 
        fileEncoding = ""
    )
    newList <- list("diffMeth" = tiles.diffMeth, "diffMeth.percent" = tiles.diffMeth.p, "gene.lookup" = df.gene.lookup)
    return(newList)
}