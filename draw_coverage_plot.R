#!/usr/bin/Rscript
# #######################################################
# #######################################################
# ### ### ### ### ### ### ### ### ### ###
# ### Draw coverage plots for Samara
# ### Date: 2020-07
# ### Aim & details:
# ### INPUT:
# ### Parametersa
# ### hard coded
# ### Output
# ### last modified
# ### TODO list/Process:
# ### ### ### ### ### ### ### ### ### ###
# #######################################################
# #######################################################

##################################
##### blabla #####
##################################
#
# samtools depth /mnt/samara/Multiplex_TSARA/WKD_jean_amaury/Run_Results/local/MaRS_test_E17MARS_jmarone_2020_07_27_1595842123/3D7-pool_S1/alignments/output_FM_SR_DD_RG.bam  > depth
# # format of "depth" tsv file: chr position depth
#
# ampliconBoudaries.tsv
# # AmpliconId primer_F primer R start ends
#
# bedfile: standard 12 column bed file
#
# R --args depth ../fasta_bed_forMars/Genes_for_samaraPipeline_2strand.bed ampliconBoudaries.tsv
# Rscript draw_coverage_plot.R depth ../fasta_bed_forMars/Genes_for_samaraPipeline_2strand.bed ampliconBoudaries.tsv


##################################
##### init envir & get input #####
##################################
tryCatch({options(width=as.numeric(strsplit(system('stty size', intern=T), ' ')[[1]][2]))}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

suppressPackageStartupMessages({
  source("./Amp_cov/plotGenes_AV_takenFromSushi.R")
})
arguments      <- commandArgs(T)

depthFile<-arguments[1]
bedFile<-arguments[2]
ampliconFile<-arguments[3]

## get & format data
###############################
depth <- read.table(file=depthFile, header=F,sep="\t", col.names=c("genes","pos","cov"), stringsAsFactors=F)
bed   <- read.table(file=bedFile,   header=F,sep="\t", col.names=c("chrom", "chromStart", "chromEnd", "genes", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"), stringsAsFactors=F)
amplicon <- read.table(file=ampliconFile,   header=T,sep="\t", stringsAsFactors=F)


plotTheGenes  <- function(currentBed, geneDinteret){
  par(mar=c(0,4,0,2)) # 5.1 4.1 4.1 2.1
  exonStarts <- as.numeric(unlist(strsplit(currentBed$blockStarts,",",fixed=T)))+1 #+1 because of plot genes
  exonEnd <- exonStarts+as.numeric(unlist(strsplit(currentBed$blockSizes,",",fixed=T)))

  genes_to_plot <- data.frame(
    chrom  = rep(currentBed$chrom,sum(currentBed$blockCount)),
    start  = exonStarts,
    stop   = exonEnd,
    gene   = rep(currentBed$genes,sum(currentBed$blockCount)),
    score  = rep(".",sum(currentBed$blockCount)),
    strand = rep(currentBed$strand,sum(currentBed$blockCount)),
    stringsAsFactors = F
  )

  # use currentBed$chromEnd+1 because the function considers that the genes goes beyond the drawing area thus do not plot. This may generate a small discrepency in scale
  pg=plotGenes(geneinfo = genes_to_plot,
               chrom = genes_to_plot$chrom[1] ,
               chromstart = currentBed$chromStart,
               chromend = currentBed$chromEnd, bheight = 0.05, bentline = FALSE, wigglefactor = 0.05, labeloffset = 0.1, fontsize = 1)
  title(ylab=paste(paste(rev(strsplit(geneDinteret,"%",fixed=T)[[1]]), collapse=" ( "),")"),font.lab=2,cex.lab=0.8,line=2.5)
  v_reduction=4*diff(par("usr")[c(3,4)])/100
  segments(x0=par("usr")[1],y0=par("usr")[3]+v_reduction,x1=par("usr")[1],y1=par("usr")[4]-v_reduction,xpd=T)
}

plotAmplicons  <- function(start, end, currentAmplicon){
  par(mar=c(2,4,0,2)) # 5.1 4.1 4.1 2.1
  plot(c(), c(), type="n", xlim=c(start, end), ylim=c(0,4)) #+1 because of plot genes

  currentAmplicon<-currentAmplicon[order(currentAmplicon$start),]
  row.names(currentAmplicon)<-NULL
  currentAmplicon$height <- (as.numeric(row.names(currentAmplicon))-1)%%3+1

  cols=c("red","blue","green");
  rect(xleft    = currentAmplicon$start,
     ybottom  = currentAmplicon$height-0.5,
     xright   = currentAmplicon$end,
     ytop    = currentAmplicon$height+0.5,
     col=cols[currentAmplicon$height], border=cols[currentAmplicon$height]
)
text(x=(currentAmplicon$start+currentAmplicon$end)/2, y = currentAmplicon$height, labels = currentAmplicon$Amplicon, cex=1 )


}

plot1GeneAmplification <- function(geneDinteret, depth, bed, amplicon){
  currentDepth    <- depth[which(depth$genes==geneDinteret),]
  currentBed      <- bed[which(bed$genes==geneDinteret),]
  currentAmplicon <- amplicon[which(amplicon$genes==geneDinteret),]

  if(nrow(currentDepth)==0){return(paste("no coverage info for",geneDinteret)) }

  if(currentBed$strand == "-"){
    newStart <- currentBed$chromEnd - currentAmplicon$end
    newEnd <- currentBed$chromEnd - currentAmplicon$start
    currentAmplicon$start <- newStart
    currentAmplicon$end <- newEnd
  }
  layout(matrix(1:3, ncol=1),height=c(1, 0.3, 0.3))
  par(mar=c(0,4,4,2)) # 5.1 4.1 4.1 2.1
  plot(currentDepth$pos, currentDepth$cov, xlim=c(currentBed$chromStart,currentBed$chromEnd), pch=20, cex=0.5)


  # gene
  plotTheGenes(currentBed, geneDinteret)

  ## amplicons
  plotAmplicons(currentBed$chromStart,currentBed$chromEnd, currentAmplicon)
}
pdf(paste(depthFile, "Kelch13.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_1343700.1%Kelch13", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"CRT.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_0709000.1%CRT", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"serine_tRNA_ligase.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_0717700.1%serine-tRNA_ligase", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"tubulin_beta_chain.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_1008700.1%tubulin_beta_chain", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"GCH1_GTP-CH.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_1224000.1%GCH1_GTP-CH", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"HRPIII.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_1372200.1%HRPIII", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"Cytb.pdf", sep = "_"))
plot1GeneAmplification("mal_mito_3%Cytb", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"DHFR-TS.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_0417200.1%DHFR-TS", depth, bed, amplicon)
graphics.off()

pdf(paste(depthFile,"PPPK-DHPS.pdf", sep = "_"))
plot1GeneAmplification("PF3D7_0810800.1%PPPK-DHPS", depth, bed, amplicon)
graphics.off()

plot3GeneAmplification <- function(vector_geneDinteret, depth, bed, amplicon){
  # code sous optimal, pas dans la logique de R, mais plus rapide à faire:
  list_Depth<-list()
  list_Bed<-list()
  list_Amplicon<-list()

  for(geneDinteret in vector_geneDinteret){
    currentDepth    <- depth[which(depth$genes==geneDinteret),]
    currentBed      <- bed[which(bed$genes==geneDinteret),]
    currentAmplicon <- amplicon[which(amplicon$genes==geneDinteret),]

    if(nrow(currentDepth)==0){return(paste("no coverage info for",geneDinteret)) }

    if(currentBed$strand == "-"){
      newStart <- currentBed$chromEnd - currentAmplicon$end
      newEnd <- currentBed$chromEnd - currentAmplicon$start
      currentAmplicon$start <- newStart
      currentAmplicon$end <- newEnd
    }
    list_Depth    <- c(list_Depth, list(currentDepth))
    list_Bed      <- c(list_Bed, list(currentBed))
    list_Amplicon <- c(list_Amplicon, list(currentAmplicon))
  }

  ylimit<-c(0,max(unlist(lapply(X=list_Depth, FUN=function(x){max(x$cov)}))))
  sizesProp<-unlist(lapply(X=list_Bed, FUN=function(x){x$chromEnd-x$chromStart}))
  layout(matrix(1:9, ncol=3,byrow=F),height=c(1, 0.2, 0.2), widths=sizesProp)
  par(mar=c(0,4,4,2)) # 5.1 4.1 4.1 2.1

  for(i in 1:3){
    par(mar=c(0,4,4,2)) # 5.1 4.1 4.1 2.1
    plot(list_Depth[[i]]$pos, list_Depth[[i]]$cov, xlim=c(list_Bed[[i]]$chromStart,list_Bed[[i]]$chromEnd), pch=20, cex=0.5, ylim=ylimit)
    # gene
    plotTheGenes(list_Bed[[i]], vector_geneDinteret[i])
    ## amplicons
    if(nrow(list_Amplicon[[i]])==0){
      plot(c(1), c(1), type="n")
    } else {
     plotAmplicons(list_Bed[[i]]$chromStart,list_Bed[[i]]$chromEnd, list_Amplicon[[i]])#+1 because of plot genes
    }
  }
}

pdf(paste(depthFile,"MDR1.pdf", sep = "_"))
plot3GeneAmplification(c("PF3D7_0522400.1%20kb_before_MDR1", "PF3D7_0523000.1%MDR1", "PF3D7_0523700.1%20kb_after_MDR1"), depth, bed, amplicon)
graphics.off()
pdf(paste(depthFile,"HRP2.pdf", sep = "_"))
plot3GeneAmplification(c("PF3D7_0831750.1%1.5kb_before_HRP2", "PF3D7_0831800.1%HRP2", "PF3D7_0832100.1%rifin_10kb_after_HRP2"), depth, bed, amplicon)
graphics.off()
pdf(paste(depthFile,"PMII.pdf", sep = "_"))
plot3GeneAmplification(c("PF3D7_1407600.1%15kb_before_PM2", "PF3D7_1408000.1%PMII", "PF3D7_1409700.1%80kb_after_PM2"), depth, bed, amplicon)
graphics.off()




# this could be used to make a better code
# https://stackoverflow.com/questions/27015449/nested-layouts-in-r
