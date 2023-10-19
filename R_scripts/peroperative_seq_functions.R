# peroperative sequencing functions

#rds file containing a colorlist 
color_translation = readRDS("~/nanocns/R-scripts/color_translation.rds")

functioncheck = function(){
  if(exists("refgenome")==F){
    print("refgenome not selected, using chm13v2.0")
    refgenome = "/home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa"}
  mega_command = paste0("/home/sturgeon/miniconda3/envs/megalodon/bin/megalodon ~/nanocns/data/functioncheck/ --outputs mods basecalls mappings --mappings-format bam ",
                        "--reference ",refgenome, " --write-mods-text --mod-motif m CG 0 --processes 10 ",
                        "--guppy-server-path /home/sturgeon/nanocns/software/guppy_v5/bin/guppy_basecall_server ",
                        "--guppy-params \"-d /home/sturgeon/nanocns/software/rerio/basecall_models/ --num_callers 2 --ipc_threads 3\"",
                        " --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --devices cuda:0 --overwrite --output-directory ~/nanocns/data/functioncheck/meg_out/ ",
                        "--suppress-progress-bars --suppress-queues-status")
  
  system("rm -r ~/nanocns/data/functioncheck/sturgeonfiles/")
  system("rm -r ~/nanocns/data/functioncheck/meg_out/")
  system("rm ~/nanocns/data/functioncheck/barcoding.tsv")
  
  system(mega_command)
  
  system(paste0("/home/sturgeon/miniconda3/bin/qcat --tsv -k RBK004 -f ",
                "~/nanocns/data/functioncheck/meg_out/basecalls.fastq > ~/nanocns/data/functioncheck//barcoding.tsv"))
  
  system("mkdir ~/nanocns/data/functioncheck/sturgeonfiles")
  system("mv ~/nanocns/data/functioncheck/meg_out/per_read_modified_base_calls.txt ~/nanocns/data/functioncheck/sturgeonfiles/per_read_modified_base_calls.txt")
  system(paste0("~/Carlo/run_sturgeon_tobed_v2.sh ~/nanocns/data/functioncheck/sturgeonfiles/ ~/nanocns/data/functioncheck/sturgeonfiles/"))
  #system2('open', args = "~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.pdf", wait = FALSE)
  res = read.table("~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.csv", sep=",", header=T)
  print(res[1:2])
  system("rm -r ~/nanocns/data/functioncheck/sturgeonfiles/")
  system("rm -r ~/nanocns/data/functioncheck/meg_out/")
  system("rm ~/nanocns/data/functioncheck/barcoding.tsv")
  print("check complete, no errors")
  
}

functioncheck_guppy = function(){
  if(exists("refgenome")==F){
    print("refgenome not selected, using chm13v2.0")
    refgenome = "/home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa"}
  
  guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i  /home/sturgeon/nanocns/data/functioncheck/ -s /home/sturgeon/nanocns/data/functioncheck/guppy_out ",
                         "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r9.4.1_e8.1_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
                         "--align_ref ", refgenome,
                         " --barcode_kits SQK-RBK004 ",
                         "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4")
  
  system(guppy_command, ignore.stdout = T, ignore.stderr=T)
  system(guppy_command)

  system(paste0("~/Carlo/sturgeon_guppy.sh ~/nanocns/data/functioncheck/guppy_out/pass/barcode02/"))
  
  #system2('open', args = "~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.pdf", wait = FALSE)
  res = read.table("~/nanocns/data/functioncheck//  ", sep=",", header=T)
  print(res[1:2])
  system("rm -r ~/nanocns/data/functioncheck/sturgeonfiles_guppy/")
  system("rm -r ~/nanocns/data/functioncheck/guppy_out/")
  print("check complete, no errors")
  
}

functioncheck_guppy_R10 = function(){
  if(exists("refgenome")==F){
    print("refgenome not selected, using chm13v2.0")
    refgenome = "/home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa"}
  
  #guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i  /home/sturgeon/nanocns/data/functioncheck_R10/ -s /home/sturgeon/nanocns/data/functioncheck/guppy_out ",
  #                       "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
  #                       "--align_ref ", refgenome,
  #                       " --barcode_kits SQK-RBK114-24 ",
  #                       "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4")
  
  guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i /home/sturgeon/nanocns/data/functioncheck_R10/ -s /home/sturgeon/nanocns/data/functioncheck_R10/guppyres ",
                         "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
                         "--disable_qscore_filtering --min_score_barcode_front 6 --min_score_barcode_rear 6 ",
                         "--align_ref ",refgenome," --barcode_kits SQK-RBK114-24 ",
                         "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  
  
  #system(guppy_command, ignore.stdout = T, ignore.stderr=T)
  system(guppy_command)
  
  system(paste0("~/Carlo/sturgeon_guppy.sh ~/nanocns/data/functioncheck_R10/guppyres/barcode02/"))
  
  #system2('open', args = "~/nanocns/data/functioncheck/sturgeonfiles/merged_probes_methyl_calls_model.pdf", wait = FALSE)
  
  system("rm -r ~/nanocns/data/functioncheck_R10/guppyres")
}

wrappert = function(main_folder,fast5, iteration, bcoverride=F){
  #this wrapper runs megalodon on a fast5 file
  #then uses qcat to split out the most frequent barcode
  #then runs sturgeon 
  #fast5 = "~/nanocns/data/example_dataII/FAR96725_b5b5d4b3_0.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  mega_command = paste0("/home/sturgeon/miniconda3/envs/megalodon/bin/megalodon ", out_folder, "  --outputs mods basecalls mappings --mappings-format bam ",
                        "--reference ",refgenome, " --write-mods-text --mod-motif m CG 0 --processes 10 ",
                        "--guppy-server-path /home/sturgeon/nanocns/software/guppy_v5/bin/guppy_basecall_server ",
                        "--guppy-params \"-d /home/sturgeon/nanocns/software/rerio/basecall_models/ --num_callers 2 --ipc_threads 3\"",
                        " --guppy-config res_dna_r941_min_modbases_5mC_CpG_v001.cfg --devices cuda:0 --overwrite --output-directory ", out_folder,"/meg_out/ ",
                        "--suppress-progress-bars --suppress-queues-status")
  
  print(mega_command)
  system(mega_command, ignore.stdout = T, ignore.stderr=T)
  
  system(paste0("mv ", out_folder,"/meg_out/* ", out_folder))
  #run qcat barcoding
  
  system(paste0("/home/sturgeon/miniconda3/bin/qcat --tsv -k RBK004 -f ",
                out_folder,"/basecalls.fastq > ",out_folder,"/barcoding.tsv"))
  #select reads from the most common barcode
  bcds = read.table(paste0(out_folder,"/barcoding.tsv"),header=T)
  if(bcoverride==F){
  
  mostfreq = names(sort(table(bcds$barcode), decreasing = T, na.last = T)[1])
  print(paste("using most frequent barcode, barcode nr",mostfreq))
  bcds = bcds[bcds$barcode==mostfreq,"name"]
  print(paste("total reads:",length(bcds)))}
  
  if(bcoverride!=F){
    print(paste("using override, barcode nr",bcoverride))
    bcds = bcds[bcds$barcode==bcoverride,"name"]
    print(paste("total reads:",length(bcds)))
  }
  
  
  #filter read calls
  methfile= read.table(paste0(out_folder,"/per_read_modified_base_calls.txt"),header=T)
  system(paste0("mv ",out_folder,"/per_read_modified_base_calls.txt ",out_folder,"/meg_out/"))
  methfile=methfile[methfile$read_id%in%bcds,]
  write.table(methfile, file=paste0(main_folder,"/per_read_modified_base_calls_it",iteration,".txt"),quote=F, row.names=F, sep="\t")
}

fix_name=function(x){
  return(gsub(gsub(x = x, pattern = "\\.\\.\\.", replacement = " - "), pattern = "\\.", replacement = " "))}

add_and_plot=function(merged_data, result_file, outfile="x"){
  
  if(exists("color_translation")==F){color_translation=readRDS("~/nanocns/R-scripts/color_translation.rds")}
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
      clr = unlist(unname(color_translation[mgd[i,"class"]]))
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

confidence_over_time_plot=function(merged_data){
  if(exists("color_translation")==F){color_translation=readRDS("~/nanocns/R-scripts/color_translation.rds")}
  mgd = merged_data
  colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
  plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="confidence over time", xlab="iteration", ylab="confidence",
       xaxt="n")
  for(i in 1:(nrow(mgd)-1)){
    clr = unlist(unname(color_translation[mgd[i,"class"]]))
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

wrappert_dorado_R10 = function(main_folder,fast5, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs dorado on a pod5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  #out_folder="~/nanocns/data/pod5_iteration/"
  
  dorado_command=paste0("~/nanocns/dorado/dorado-0.3.2-linux-x64/bin/dorado basecaller ~/nanocns/models/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ", out_folder,
                       " --reference /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --modified-bases 5mCG_5hmCG > ", out_folder,"/mapped.bam")
  print("running dorado")
  system(dorado_command)
  
  system(paste0("samtools sort -@ 10 ", out_folder, "/mapped.bam -o ",out_folder, "/mapped_srt.bam"))
  system(paste0("rm ",out_folder, "/mapped.bam"))
  
  modkit_adjust = paste0("~/nanocns/modkit/modkit adjust-mods --convert h m ", out_folder, "/mapped_srt.bam ",out_folder, "/mapped_srt_adj.bam")
  system(modkit_adjust)
  system(paste0("rm ",out_folder, "/mapped_srt.bam"))
  system(paste0('samtools index ', out_folder, "/mapped_srt_adj.bam"))
  print("modkit extract")
  modkit_extract = paste0("~/nanocns/modkit/modkit extract ",out_folder, "/mapped_srt_adj.bam ",out_folder,"/modkit_extract.txt" )
  system(modkit_extract)
  print("running sturgeon")
  
  system(paste0("~/Carlo/run_sturgeon_tobed_modkit.sh ",out_folder, " ",out_folder))
  
 # barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
 # system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))
 # 
 # #include unclassified reads yes or no?
 # if(include_unclassified==T){
 #   system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_it",iteration,"_unclassified.bam"))}

}

wrappert_guppy_R10 = function(main_folder,fast5, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  
 #guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i ",out_folder, " -s ", out_folder,"/guppy_out ",
 #                       "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
 #                       "--disable_qscore_filtering --min_score_barcode_front 6 --min_score_barcode_rear 6 ",
 #                       "--align_ref /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --barcode_kits SQK-RBK114-24 ",
 #                       "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
 
 guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6.4.6/ont-guppy_6.4.6_linux64/ont-guppy/bin/guppy_basecaller -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                        "-c /home/sturgeon/nanocns/software/guppy_v6.4.6/ont-guppy_6.4.6_linux64/ont-guppy/data/dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
                        "--disable_qscore_filtering --min_score_barcode_front 6 ",
                        "--align_ref /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --barcode_kits SQK-RBK114-24 ",
                        "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
 
  #guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i ",out_folder, " -s ", out_folder,"/guppy_out ",
  #                       "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r10.4.1_e8.2_400bps_modbases_5mc_cg_fast.cfg --device auto --bam_out ",
  #                       "--disable_qscore_filtering --min_score_barcode_front 6 --min_score_barcode_rear 6 ",
  #                       "--align_ref /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --barcode_kits SQK-RBK114-24 ",
  #                       "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4")
  #system(guppy_command, ignore.stdout = T, ignore.stderr=T)
  system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))
  
  #include unclassified reads yes or no?
  if(include_unclassified==T){
  system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_it",iteration,"_unclassified.bam"))}
  
  system(paste0("~/Carlo/sturgeon_guppy.sh ",main_folder))
  
  
}

#needs to be tested still
wrappert_guppy_R10_guppy6.5 = function(main_folder,fast5, iteration, bcoverride=F, include_unclassified=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  

  guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6.5.7/ont-guppy/bin/guppy_basecaller -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                         "-c ~/nanocns/software/guppy_v6.5.7/ont-guppy/data/dna_r10.4.1_e8.2_400bps_5khz_modbases_5mc_cg_hac_mk1c.cfg --device auto --bam_out ",
                         "--disable_qscore_filtering --min_score_barcode_front 6 ",
                         "--align_ref /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --barcode_kits SQK-RBK114-24 ",
                         "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4 --allow_inferior_barcodes")
  
  system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  system(paste0("mv ", out_folder,"/guppy_out/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))
  
  #include unclassified reads yes or no?
  if(include_unclassified==T){
    system(paste0("mv ", out_folder,"/guppy_out/unclassified/*.bam ",main_folder,"/guppy_output_it",iteration,"_unclassified.bam"))}
  
  system(paste0("~/Carlo/sturgeon_guppy.sh ",main_folder))
  
  
}

wrappert_guppy_R9 = function(main_folder,fast5, iteration, bcoverride=F){
  #this wrapper runs guppy on a fast5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  out_folder = paste0(main_folder,"/iteration_",iteration)
  
  system(paste0("mkdir ",out_folder))
  system(paste0("cp ", fast5," ", out_folder))
  
  
  guppy_command = paste0("/home/sturgeon/nanocns/software/guppy_v6/bin/guppy_basecaller -i ",out_folder, " -s ", out_folder,"/guppy_out ",
                         "-c /home/sturgeon/nanocns/software/guppy_v6/data/dna_r9.4.1_e8.1_modbases_5mc_cg_hac.cfg --device auto --bam_out ",
                         " --disable_qscore_filtering ",
                         "--align_ref /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --barcode_kits SQK-RBK004 ",
                         "--chunk_size 2000 --chunks_per_runner 256 --chunks_per_caller 10000 --gpu_runners_per_device 4  --num_base_mod_threads 4")
  
  system(guppy_command, ignore.stdout = T, ignore.stderr=T)
  #system(guppy_command)
  
  barcode = ifelse(nchar(bcoverride)==1, paste0("0", bcoverride), bcoverride)
  system(paste0("mv ", out_folder,"/guppy_out/pass/barcode",barcode,"/*.bam ",main_folder,"/guppy_output_it",iteration,".bam"))
  system(paste0("~/Carlo/sturgeon_guppy.sh ",main_folder))

  
}

library(QDNAseq)
library(HMM)

bins=readRDS("~/nanocns/chm13_v2.0_bins.rds")

plot_cnv_from_bam = function(bam){
  readCounts <- binReadCounts(bins, bamfiles=bam)
  exportBins(readCounts, file="~/nanocns/bins_sample_counts.bed", format="bed", filter=F, logTransform=F)
  reference=read.table("~/nanocns/bins_export.bed", skip = 1)
  colnames(reference) = c("chrom", "start", "end", "name", "coverage", "orientation")
  sample = read.table("~/nanocns/bins_sample_counts.bed", skip = 1)
  colnames(sample) = c("chrom", "start", "end", "name", "coverage", "orientation")
  system("rm ~/nanocns/bins_sample_counts.bed")
  relevant_genes = read.table("~/nanocns/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed", header=T)

  totreads = sum(sample$coverage)
  sample$refcov = reference$coverage
  sample = sample[sample$chrom!="Y"&sample$chrom!="X",]
  sample$reltoref = sample$coverage/sample$refcov
  sample=na.omit(sample)
  blacklistbins = sample[sample$refcov<500,]
  blacklistbins = row.names(blacklistbins)
  sample = sample[!row.names(sample)%in%blacklistbins,]
  mn = mean(sample$reltoref)
  totreads=sum(sample$coverage)
  sample$logtr = log(sample$reltoref/mn, 2)
  maxy= ifelse(max(sample$logtr)>3, max(sample$logtr),3)

  plot("",xlim=c(0,nrow(sample)), ylim=c(-3,maxy), xlab="position",xaxt="n", ylab="log ratio", main=paste0("CNV plot, nreads=",totreads))

  ##points(x=1:nrow(sample), y=sample$logtr, pch=16, cex=0.5)

  #abline(h=0, col="red")

  for(i in 1:(nrow(sample)-1)){if(sample[i,"chrom"]!=sample[i+1, "chrom"]){
    abline(v=i)
    text(x = i, adj=c(1,0), labels = sample[i,"chrom"], y = maxy)
  }}



  sample$categ = round(sample$logtr, 0)
  sample$categ = ifelse(sample$categ<(-1), -1, ifelse(sample$categ >1, 1, sample$categ))
  emiprobs = matrix(c(.73,.25,.02, 0.25,.4,.25,.02,.25,.73),3)

  hmm = initHMM(c("-1","0","1"), c("-1","0","1"), transProbs=matrix(c(.98,.01,.01,0.01,0.98,0.01, 0.01,0.01,0.98),3), emissionProbs=emiprobs)

  res=data.frame()
  for(chr in unique(sample$chrom)){
    chrres = HMM::posterior(hmm, as.character(sample[sample$chrom==chr,"categ"]))
    fres = t(chrres)

    res=rbind(res, fres)
  }

 dfres=data.frame(res)

  row.names(dfres)=1:nrow(dfres)
  dfres$max_col <- apply(dfres, 1, which.max)
  dfres$plot= ifelse(dfres$max_col ==1 , -1, ifelse(dfres$max_col==2, 0, 1))
  dfres$plot = ifelse(dfres[,1]>0.9, -1,
                      ifelse(dfres[,2]>0.9, 0,
                             ifelse(dfres[,3]>0.9, 1, 0)))
  dfres$color= ifelse(dfres[,1]>0.9, "blue",
                      ifelse(dfres[,2]>0.9, "darkgrey",
                             ifelse(dfres[,3]>0.9, "green","grey")))

  #plot("",xlim=c(0,nrow(sample)), ylim=c(-3,max(sample$logtr)))
  for(i in 1:nrow(relevant_genes)){
    dfres[relevant_genes[i, "binnr"],"color"]="red"
  }


  dfres[blacklistbins,"color"]="white"
  points(x=1:nrow(sample), y=sample$logtr, pch=16, cex=0.5, col=dfres$color)

  lines(x=1:nrow(dfres), y=dfres$plot, pch=16, col="black", lwd=2)

  for(i in 1:nrow(relevant_genes)){
    ycoord = sample[relevant_genes[i,"binnr"],"logtr"]
    color = ifelse(dfres[relevant_genes[i,"binnr"],"color"])
    text(x=relevant_genes[i,"binnr"]+3, y=ycoord-0.1+relevant_genes[i,"yoffset"], cex=0.5, srt=-90,adj=c(0.5,1), labels = relevant_genes[i, "name"], col=color, pos=1)
  }

}
plot_cnv_from_live_dir=function(directory){
   #directory = "~/nanocns/data/pmc_live_12/"
   samtools = "~/nanocns/software/samtools-1.17/samtools"
   
   if(dir.exists(paste0(directory,"/merged_bams/"))==F){system(paste0("mkdir ", directory,"merged_bams/"))}else{
     system(paste0("rm ", directory,"/merged_bams/merged_bam.bam"))
   }
   
   system(paste0(samtools, " merge -@ 10 -O BAM -o ",directory,"/merged_bams/merged_bam.bam ", directory,"/*.bam" ))
   
   plot_cnv_from_bam(paste0(directory,"/merged_bams/merged_bam.bam"))
   
 }

library(depmixS4)
plot_cnv_from_bam_depmix=function(bam){
  readCounts <- binReadCounts(bins, bamfiles=bam, isNotPassingQualityControls=NA)
  
  exportBins(readCounts, file="~/nanocns/bins_sample_counts_tmp.bed", format="bed", filter=F, logTransform=F)
  
  reference=read.table("~/nanocns/bins_export.bed", skip = 1)
  colnames(reference) = c("chrom", "start", "end", "name", "coverage", "orientation")
  sample = read.table("~/nanocns/bins_sample_counts_tmp.bed", skip = 1)
  system("rm ~/nanocns/bins_sample_counts_tmp.bed")
  colnames(sample) = c("chrom", "start", "end", "name", "coverage", "orientation")
  
  totreads = sum(sample$coverage)
  sample$refcov = reference$coverage
  sample = sample[sample$chrom!="Y"&sample$chrom!="X",]
  sample$reltoref = sample$coverage/sample$refcov
  sample=na.omit(sample)
  blacklistbins = sample[sample$refcov<500,]
  blacklistbins = row.names(blacklistbins)
  sample = sample[!row.names(sample)%in%blacklistbins,]
  mn = mean(sample$reltoref)
  totreads=sum(sample$coverage)
  sample$logtr = log(sample$reltoref/mn, 2)
  maxy= ifelse(max(sample$logtr)>3, max(sample$logtr),3)
  
  trprobs = rep(0.01, times=25)
  trprobs[1] = 0.96
  last=1
  for(state in 2:5){
    trprobs[last+6]=0.96
    last=last+6
  }
  
  
  sample=sample[sample$coverage!=0,]
  mod <- depmix(response = logtr ~ 1, data = sample, nstates = 5, trstart = trprobs,respstart=c(-1,0.3, -0.5, 0.3, 0, 0.3, 0.5, 0.3, 1, 0.3))
  #fm <- fit(mod)
  fm <- mod
  summary(fm)
  statemeans= data.frame(summary(fm))
  statemeans=statemeans$Re1..Intercept.
  esttrans <- posterior(fm, type="viterbi" )
  transtable=data.frame(states=length(statemeans), means=statemeans)
  translate = function(x){return(transtable[x,"means"])}
  
  esttrans$corrmn = sapply(esttrans$state,translate)
  sample$hmmpred=esttrans$corrmn
  sample$color = ifelse(sample$hmmpred<(-0.25),"blue",
                        ifelse(sample$hmmpred>0.25,"green", "grey"))
  relevant_genes = read.table("~/nanocns/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed", header=T)
  
  #sample[relevant_genes$binnr,"color"]="red"
  
  plot("",xlim=c(0,max(as.numeric(row.names(sample)))), ylim=c(-3,maxy), xlab="position",xaxt="n", ylab="log ratio", main=paste0("CNV plot, nreads=",totreads))
  
  
  points(x=as.numeric(row.names(sample)),y=sample$logtr, type="p", pch=16, col=sample$color)
  #lines(x = as.numeric(row.names(sample)), y=esttrans$corrmn, col="red" )
  
  
  segments(x0 = as.numeric(row.names(sample))-0.5, x1=as.numeric(row.names(sample))+0.5, y0=esttrans$corrmn,y1=esttrans$corrmn, col="red", lwd=2 )
  
  for(i in 1:(nrow(sample)-1)){if(sample[i,"chrom"]!=sample[i+1, "chrom"]){
    abline(v=row.names(sample)[i])
    text(x = row.names(sample)[i], adj=c(1,0), labels = sample[i,"chrom"], y = maxy)
  }}
  
  for(i in 1:nrow(relevant_genes)){
    chrnr = strsplit(relevant_genes[i,"chrom"], "chr")[[1]][2]
    gstart = relevant_genes$start[i]
    ycoord = sample[sample$chrom==chrnr&sample$start<gstart&sample$end>gstart,"logtr"]
    xcoord = row.names(sample[sample$chrom==chrnr&sample$start<relevant_genes$start[i]&sample$end>relevant_genes$start[i],])
    xcoord=as.numeric(xcoord)
    points(x=xcoord, y=ycoord, col="red")
    text(x=xcoord+3, y=ycoord-0.1+relevant_genes[i,"yoffset"], cex=0.5, srt=-90,adj=c(0.5,1), labels = relevant_genes[i, "name"], pos=1)
  }
}
plot_cnv_from_live_dir_depmix=function(directory){
  #directory = "~/nanocns/data/pmc_live_12/"
  samtools = "~/nanocns/software/samtools-1.17/samtools"
  
  if(dir.exists(paste0(directory,"/merged_bams/"))==F){system(paste0("mkdir ", directory,"merged_bams/"))}else{
    system(paste0("rm ", directory,"/merged_bams/merged_bam.bam"))
  }
  
  system(paste0(samtools, " merge -@ 10 -O BAM -o ",directory,"/merged_bams/merged_bam.bam ", directory,"/*.bam" ))
  
  plot_cnv_from_bam_depmix(paste0(directory,"/merged_bams/merged_bam.bam"))
  
}

library(DNAcopy)

plot_cnv_from_bam_DNAcopy = function(bam, makeplot=T, lines_only=F, binsize=1e6){
  
  
  readCounts <- binReadCounts(bins, bamfiles=bam, isNotPassingQualityControls=NA, minMapq=2, isDuplicate=NA, isSecondaryAlignment=NA)
  exportBins(readCounts, file="~/nanocns/bins_sample_counts_tmp.bed", format="bed", filter=F, logTransform=F)
  reference=read.table("~/nanocns/bins_export.bed", skip = 1)
  colnames(reference) = c("chrom", "start", "end", "name", "coverage", "orientation")
  sample = read.table("~/nanocns/bins_sample_counts_tmp.bed", skip = 1)
  system("rm ~/nanocns/bins_sample_counts_tmp.bed")
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
      relevant_genes = read.table("~/nanocns/relevant_genes_with_chm13v2_1Mb_bin_nrs.bed", header=T)
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
plot_cnv_from_live_dir_DNAcopy=function(directory){
  #directory = "~/nanocns/data/pmc_live_12/"
  samtools = "~/nanocns/software/samtools-1.17/samtools"
  
  if(dir.exists(paste0(directory,"/merged_bams/"))==F){system(paste0("mkdir ", directory,"/merged_bams/"))}else{
    system(paste0("rm ", directory,"/merged_bams/merged_bam.bam"))
  }
  
  system(paste0(samtools, " merge -@ 10 -O BAM -o ",directory,"/merged_bams/merged_bam.bam ", directory,"/*.bam" ))
  
  out=plot_cnv_from_bam_DNAcopy(paste0(directory,"/merged_bams/merged_bam.bam"))
  
}


dorado_CNV_plotting = function(pod5_folder, out_folder){
  #this wrapper runs dorado on a pod5 file
  #then runs sturgeon on a selected barcode
  #iteration = 3
  #bcoverride = 1
  #main_folder = "~/nanocns/data/guppytest"
  #system(paste0("mkdir ", main_folder))
  #fast5 = "/var/lib/minknow/data/AMC_run_1/no_sample/20230320_1146_MN40017_FAV81912_63a68867/fast5/FAV81912_63a68867_c5a6d74a_2.fast5"
  
  system(paste0("mkdir ",out_folder))
  
  #out_folder="~/nanocns/data/pod5_iteration/"
  
  dorado_command=paste0("~/nanocns/dorado/dorado-0.3.2-linux-x64/bin/dorado basecaller ~/nanocns/models/dna_r10.4.1_e8.2_400bps_hac@v4.2.0 ", pod5_folder,
                        " --reference /home/sturgeon/nanocns/data/chm13_V2/chm13v2.0.fa --modified-bases 5mCG_5hmCG > ", out_folder,"/mapped.bam")
  print("running dorado")
  system(dorado_command)
  
  system(paste0("samtools sort -@ 10 ", out_folder, "/mapped.bam -o ",out_folder, "/mapped_srt.bam"))
  system(paste0("samtools index ", out_folder, "/mapped_srt.bam"))
  system(paste0("rm ",out_folder, "/mapped.bam"))
  out=plot_cnv_from_bam_DNAcopy(paste0(out_folder,"/mapped_srt.bam"))
  
}

