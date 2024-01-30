# peroperative sequencing functions

library(QDNAseq)
library(DNAcopy)
#required files

##windows format paths
#sturgeon rscripts dir, windows format
rscriptsdir = "c:/Users/sturgeon/Desktop/sturgeon/R_scripts/"

#rds file containing a colorlist 

color_translation = readRDS(paste0(rscriptsdir, "/color_translation.rds"))
#this is a temporary file used for plotting the cnv plots, it will be generated and deleted, 
#you only have to make sure the location is accessible
cnvplot_tmpfile = paste0(rscriptsdir, "/bins_sample_counts_tmp.bed")
#this is the reference file used for plotting the cnv plots
cnvplot_reffile = paste0(rscriptsdir, "/bins_export.bed")

#DNAcopy bins
bins=readRDS(paste0(rscriptsdir, "/chm13_v2.0_bins.rds"))
relevant_genes_cnvplot = paste0(rscriptsdir, "/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed")

#linux formats
#point to the samtools location as reconised within the wsl (so leading with /mnt/)
samtools_location = "/mnt/c/Users/sturgeon/Desktop/sturgeon_files/samtools/samtools-1.19.2/samtools" 
#make sure the access rights are set for this sh script
sturgeon_shell = "/mnt/c/Users/sturgeon/Desktop/sturgeon/R_scripts/run_sturgeon_windows.sh"                                


#below this everything should work without changes

#fix weird errors popping up from the color translation
color_translation_frame = data.frame(names = names(color_translation), color_translation)
row.names(color_translation_frame)=1:nrow(color_translation_frame)


fix_name=function(x){
  return(gsub(gsub(x = x, pattern = "\\.\\.\\.", replacement = " - "), pattern = "\\.", replacement = " "))}

add_and_plot_windows=function(merged_data, result_file, outfile="x"){
  
  if(exists("color_translation")==F){color_translation=readRDS("R-scripts/color_translation.rds")}
  reslist_2 = t(data.frame(read.table(result_file, sep=",")))
  reslist_2=reslist_2[2:nrow(reslist_2),]
  colnames(reslist_2)=c("class", "score")
  tm = Sys.time()
  tm= paste(format(as.POSIXct(tm), format = "%H:%M"))
  reslist_2=rbind(reslist_2,data.frame(class="TIME", score=tm))
  
  if(nrow(merged_data)==0){
    return(reslist_2)
  }else{
    
    mgd = merge(merged_data, reslist_2, by="class") 
    colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
    
    plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
         xaxt="n") 
    for(i in 1:(nrow(mgd)-1)){
      clr = color_translation_frame[color_translation_frame$names==mgd[i,"class"],2]
      lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
      if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
    }
    maxtics= ifelse(ncol(mgd)>5,ncol(mgd), ncol(mgd)+5)
    axis(side = 1, at = 0:maxtics, labels = 1:(maxtics+1))
    axis(side = 1, at = 0:(ncol(mgd)-2), labels = mgd[mgd$class=="TIME",2:ncol(mgd)], line=1, tick=F, cex.axis=0.7)
    abline(h=0.95, col="red")
    abline(h=0.80, col="orange")
    return(mgd)
  }
}

confidence_over_time_plot_windows=function(merged_data, color_translation= color_translation_frame){
  mgd = merged_data
  colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
  plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
       xaxt="n")
  for(i in 1:(nrow(mgd)-1)){
    clr = color_translation_frame[color_translation_frame$names==mgd[i,"class"],2]
    lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
    if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
  }
  maxtics= ifelse(ncol(mgd)>5,ncol(mgd), ncol(mgd)+5)
  axis(side = 1, at = 0:maxtics, labels = 1:(maxtics+1))
  axis(side = 1, at = 0:(ncol(mgd)-2), labels = mgd[mgd$class=="TIME",2:ncol(mgd)], line=1, tick=F, cex.axis=0.7)
  abline(h=0.95, col="red")
  abline(h=0.80, col="orange")
  #return(mgd)
}

plot_cnv_from_bam_DNAcopy = function(bam, makeplot=T, lines_only=F, binsize=1e6){
  readCounts <- binReadCounts(bins, bamfiles=bam, isNotPassingQualityControls=NA, minMapq=2, isDuplicate=NA, isSecondaryAlignment=NA)
  exportBins(readCounts, file=cnvplot_tmpfile, format="bed", filter=F, logTransform=F)
  reference=read.table(cnvplot_reffile, skip = 1)
  colnames(reference) = c("chrom", "start", "end", "name", "coverage", "orientation")
  sample = read.table(cnvplot_tmpfile, skip = 1)
  system(paste0("rm ",cnvplot_tmpfile))
  colnames(sample) = c("chrom", "start", "end", "name", "coverage", "orientation")
  totreads_title = sum(sample$coverage)
  #print(totreads_title)
  
  sample$coverage = as.numeric(sample$coverage)+0.001
  totreads = sum(sample$coverage)
  sample$refcov = reference$coverage
  sample = sample[sample$chrom!="Y"&sample$chrom!="X",]
  sample$reltoref = sample$coverage/sample$refcov
  ##sample=na.omit(sample)
  blacklistbins = sample[sample$refcov<500,]
  blacklistbins = row.names(blacklistbins)
  sample = sample[!row.names(sample)%in%blacklistbins,]
  #sample=sample[sample$coverage>0,]
  mn = mean(sample$reltoref)
  totreads=sum(sample$coverage)
  sample$logtr = log(sample$reltoref/mn, 2)
  sample$logtr=ifelse(sample$logtr>3,3, ifelse(sample$logtr<(-3),-3, sample$logtr))
  sample$chromname=ifelse(nchar(sample$chrom)==1, paste0("chr0", sample$chrom), paste0("chr", sample$chrom))
  bam_cna = CNA(sample$logtr, chrom=sample$chromname, maploc=sample$start+binsize, data.type="logratio", sampleid = "cnvplot")
  smoothed.CNA.object <- smooth.CNA(bam_cna)
  segment.smoothed.CNA.object <- segment(smoothed.CNA.object, verbose=1)
  #plot(segment.smoothed.CNA.object, plot.type="w")
  pointcolorframe = data.frame(pointx=1:length(bam_cna$cnvplot), pointy=bam_cna$cnvplot, pointcolor="black")
  for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]>0.5){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="green"
    }
    if(segment.smoothed.CNA.object$output[x, "seg.mean"]<(-0.5)){
      pointcolorframe[segment.smoothed.CNA.object$segRows[x,1]:segment.smoothed.CNA.object$segRows[x,2],"pointcolor"]="blue"
    }    }
  if(makeplot){
    plot("", xlim=c(0,length(bam_cna$cnvplot)), ylim=c(-3,3), main=paste0("CNVs nreads=",totreads_title), xaxt="n", xlab="genomic pos",
         ylab="log ratio")
    if(lines_only==F){
      points(x=pointcolorframe$pointx, y=pointcolorframe$pointy, pch=16, cex=0.5, col=pointcolorframe$pointcolor)
    }
    for(x in 1:nrow(segment.smoothed.CNA.object$segRows)){
      segments(x0=  segment.smoothed.CNA.object$segRows[x,1], x1=segment.smoothed.CNA.object$segRows[x,2],
               y0=segment.smoothed.CNA.object$output[x,"seg.mean"], y1=segment.smoothed.CNA.object$output[x,"seg.mean"],
               col="red", lwd=2)
    }
    fr = data.frame(chrom=smoothed.CNA.object$chrom, loc=smoothed.CNA.object$maploc, val=smoothed.CNA.object$cnvplot)
    nr = 1
    for(i in 1:(nrow(fr)-1)){
      if(fr$chrom[i]!=fr$chrom[i+1]){
        abline(v=as.numeric(row.names(fr)[i])+0.5, col="grey")
        chrname= gsub(pattern = "chr", replacement = "", x = fr[i,"chrom"])
        ypos = ifelse(nr %% 2 == 1, 2, 1.8)
        nr=nr+1
        text(x = i, adj=c(1,0), labels = chrname, y = ypos)
      }
    }
    if(lines_only==F){
      relevant_genes = read.table(relevant_genes_cnvplot, header=T)
      fr$start=fr$loc-binsize
      for(i in 1:nrow(relevant_genes)){
        genestart = relevant_genes[i,"start"]
        genechrom = relevant_genes[i,"chrom"]
        genechrom = strsplit(genechrom, split="hr")[[1]][2]
        genechrom=ifelse(nchar(genechrom)==1, paste0("chr0", genechrom), paste0("chr", genechrom))
        binline = fr[fr$chrom==genechrom&fr$start<genestart&(fr$start+binsize*2)>genestart,]
        ycoord = binline$val
        text(x=as.numeric(row.names(binline)), y=ycoord-0.1+relevant_genes[i,"yoffset"], cex=0.5, srt=-90,adj=c(0.5,1), labels = relevant_genes[i, "name"], pos=1)
        points(x=as.numeric(row.names(binline)), y=ycoord, col="red")
      }}
  }
  output = list(pointdata=pointcolorframe, segdata=segment.smoothed.CNA.object)
  return(output)
}



