loadWorkspace <- function(){

    list.of.packages <- c("data.table", "dplyr", "plyr", "circlize", "tibble")

    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if(length(new.packages)>0) {
        install.packages(new.packages)
    }

    lapply(list.of.packages, require, character.only = TRUE)

    list.of.packages <- c("genomation", "methylKit", "GenomicRanges", "ggbio", "RColorBrewer", "BSgenome.Mmusculus.UCSC.mm10.masked", "ComplexHeatmap", "compEpiTools")

    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

    if(length(new.packages)>0) {
        source("http://www.bioconductor.org/biocLite.R")
        biocLite(new.packages)
    }

    lapply(list.of.packages, require, character.only = TRUE)

}