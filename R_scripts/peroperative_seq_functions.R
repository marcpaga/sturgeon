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

add_and_plot=function(merged_data, result_file){
  
  if(exists("color_translation")==F){color_translation=readRDS("~/nanocns/R-scripts/color_translation.rds")}
  reslist_2 = t(data.frame(read.table(result_file, sep=",")))
  reslist_2=reslist_2[2:nrow(reslist_2),]
  colnames(reslist_2)=c("class", "score")
  
  if(nrow(merged_data)==0){
    return(reslist_2)
  }else{
    
    mgd = merge(merged_data, reslist_2, by="class") 
    colnames(mgd)[ncol(mgd)]=paste0("iteration_", ncol(mgd)-1)
    plot('', xlim=c(0,ifelse(ncol(mgd)>5, ncol(mgd)+5, 10)), ylim=c(0,1), main="performance over time", xlab="iteration", ylab="confidence") 
    for(i in 1:nrow(mgd)){
      clr = unlist(unname(color_translation[mgd[i,"class"]]))
      lines(x = 0:(ncol(mgd)-2), y= as.numeric(mgd[i,2:ncol(mgd)]), col=clr)
      if(as.numeric(mgd[i,ncol(mgd)])>0.5){text(x=ncol(mgd), y=as.numeric(mgd[i,ncol(mgd)])-0.02, labels = mgd[i,"class"])}
    }
    abline(h=0.95, col="red")
    abline(h=0.80, col="orange")
    return(mgd)}
}
