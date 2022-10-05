if(FALSE){
  closeAllConnections()
  assign("last.warning", NULL, envir = baseenv())
  rm(list = ls())
  closeAllConnections()
  setwd("~/github/procSC/example/flames")
}
{  
library(lmerTest)
library(jsonlite)
library(DirichletReg)
  currdir = getwd()
  
params_file="params.json"
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) params_file = args[1]


print(params_file)
params = read_json(params_file); options(params)  
if(!is.null(getOption("include"))){
  params1 = read_json(getOption("include")); options(params1); options(params);
}
try(source(paste(getOption("fspls_home","~/bitbucket/fs-pls/"),"R/projection_free_functions.R",sep="/")))
try(source(paste(getOption("procSC_home","~/procSC"),"R/procSC_lib.R",sep="/")))


outdir = getOption("outdir","results")
if(length(args)>1) {
 outdir = args[2]
}else{
 outdir=paste0(outdir,format(Sys.time(), "%a_%b_%d_%H_%M_%S_%Y"))
}
options("outdir"=outdir)
working_dir = getOption("working_dir","./")



#try(source("/data/gpfs/projects/punim1466/sw/multigo/R/procSC_lib.R"))

whitelist = NULL
if(!is.null(getOption("transcript_whitelist"))){
  whitelist = read.csv(getOption("transcript_whitelist"),head=F)[,1]
}
defs_=.readDefs(paste(currdir,getOption("defs","defs.json"),sep="/"))

setwd(working_dir)
samples=read.csv(getOption("meta_file","samples.txt")[[1]],sep="\t",header=T)
##NOW INTEGRATE ANY CLINICAL DATA WITH SAMPLE TABLE
#clinical_data = getOption("clinical_data",NULL)
#if(!is.null(clinical_data)){
#  samples=.addClinicalData(samples, read.table(clinical_data[[1]], sep="\t", header=T))
#}
##FOLLOWING IF WE DONT HAVE TOTAL READS COUNTS AS COLUMN IN meta_file
if(nchar(getOption("sumsInMeta",""))==0 || is.null(samples[getOption("sumsInMeta")])){
  samples=.addCountData(samples)
}
samples=.addSplitColumns(samples, getOption("cohortSplitName"),getOption("cohortSplit"))
samples=.addMergeColumns(samples, getOption("mergeColumns"),getOption("cohortSplit"))

samples=.addSumColumns(samples, "pseudobulk")
samples = .addSumColumns(samples, "library")

##REMOVING CELLS WITH LOW COUNT PER CELL
samples =samples[ samples[[getOption("sumsInMeta")]]>=getOption("min_reads_per_cell" ,1000),,drop=F]


connection = gzfile( getOption("input_file"),"rb")
header=strsplit(readLines(connection,1),getOption("split",",")[[1]])[[1]]


samples1 = lapply(defs_, function(def){
  print(attr(def,"nme"))
  mi = match( names(def$case),names(samples))
  names(mi)=names(def$case)
  if(length(which(is.na(mi)))>0) stop("paste problem with defs ")
  samp1 = .getDataFrame(def, samples, header)
  adj_random = getOption("adjust_random",NULL)
  ##this adds a replicate column ,which is just the index of the adj_random factor
  if(is.null(samp1$Replicate) && length(adj_random)>0){
    inds_t = lapply(unique(samp1$casecontrol), function(t) which(samp1$casecontrol==t))
    replicate =data.frame( matrix(NA,nrow = nrow(samp1), ncol = length(adj_random)))
    for(j in 1:length(adj_random)){
       for(ss in inds_t){
         replicate[ss,j] = as.numeric(factor(samp1[[adj_random[[j]]]][ss]))
       }
    }
    names(replicate)=paste(unlist(adj_random),"replicate",sep=".")
    samp1 = cbind(samp1, replicate)
  }
  
  if(nrow(samp1)==0) stop(" somthing wrong with definition, no barcodes selected")
  tab = table(samp1$casecontrol)
  print(tab)
  samp2 = .splitFrameMultiple(def, samp1)
  min_cell_count=unlist(lapply(samp2, function(x) {
    t=table(x$casecontrol); 
    if(length(t)==1) 0 else min(t) ;
  })
    ,rec=F)
  samp2 = samp2[min_cell_count>=getOption("min_cells",10)]
  if(length(samp2)==0) return(NULL)
  psB=lapply(samp2,.pseudoBulk, def)
  psB = psB[ !unlist(lapply(psB, is.null))]   
  psB
})
samples1 =samples1[ unlist(lapply(samples1, length))>0]   
#df_all
#"NAME:Cell_Type:Cell_State:Cohort:Patient:Count:Sex:Variant:Time:Age:Cohort_State:Patient_Cohort:ReadCount_pseudobulk:CellCount_pseudobulk:AvgRead_pseudobulk:ReadCount_library:CellCount_library:AvgRead_library:index:index_header:casecontrol:all:xT:xG:index_within"
#df_summ
#"Cell_Type:Cell_State:Cohort:Patient:Count:Sex:Variant:Time:Age:Cohort_State:Patient_Cohort:ReadCount_pseudobulk:CellCount_pseudobulk:AvgRead_pseudobulk:ReadCount_library:CellCount_library:AvgRead_library:casecontrol:all:xT:xG:FALSE"
if(length(samples1)==0) stop(" problem. not enough cells")

#psb1=samples1$VIC_infected_vs_uninfected_child$all
dir.create(outdir)

outf=vector("list", length(samples1))
header_line=
  getOption("header",
  "ID:Name:Cell:pv_cell_fixed:pv_cell_rand:wilcox_pv_cell:dir_pv_cell:pv_pb_fixed:pv_pb_rand:wilcox_pv_pb:dir_pv_pb:beta_cell_fixed:beta_cell_rand:dir_beta_cell:beta_pb_fixed:beta_pb_rand:dir_beta_pb"
  )
#  dir_pv_pb dir_beta_pb
header1=strsplit(header_line, ":")[[1]]
                
for(j in 1:length(samples1)){
  outf[[j]] = gzfile(paste(outdir, "/",names(samples1)[[j]],".tsv.gz",sep=""),open="w")
  writeLines(paste(header1,collapse="\t" ), outf[[j]])
}

nextLine=strsplit(readLines(connection,1),getOption("split")[[1]])[[1]]

cnt=0; cnt1 = 0;
max_count=getOption("max_cnt",1e9)
min_num_reads = getOption("min_num_reads")

transcript_column = getOption("transcriptColumn",NULL)

target_genes = unlist(getOption("target_genes",c()))
gene_base_mean_thresh = getOption("gene_base_mean_thresh",0)
transcript_base_mean_thresh = getOption("transcript_base_mean_thresh",0)
min_num_transcripts = getOption("min_num_transcripts",2)
min_total_sum = getOption("min_total_sum",200)
var_thresh = getOption("var_thresh",0.001)
while(cnt<max_count){
  print(paste(cnt, cnt1,nextLine[[1]]))
  
  matr=.readMatrix(connection,nextLine)
  cnt=cnt+1
  cnt1 = cnt1+length(matr)
  nextLine = attr(matr,"next")
  if(length(target_genes)>0){
   if(!(attr(matr,"gene") %in% target_genes)){
     next;   
    }else{
      #stop("!!")
    }
  }
 
  to_include =rep(T, length(matr))
  if(!is.null(transcript_column) && !is.null(whitelist)){
    transcript_ids = unlist(lapply(names(matr), function(m1){
      strsplit(m1,"\\|")[[1]][transcript_column]
    }))
    to_include = (!is.na(match(transcript_ids, whitelist)))
  }
   if(length(which(to_include))<min_num_transcripts) next;
  matr_include = matr[to_include]
  gene = attr(matr,"gene")
  pv_res_all =  lapply(samples1,function(psB){
      results1 = lapply(psB, function(psb1){
        matr1 = data.frame(lapply(matr_include, function (line) as.numeric(line[psb1$df_all$index_header])))
        matr2 =  t(data.frame(lapply(psb1$sep, function(sep1){
          apply(matr1[sep1$index_within,,drop=F],2,sum)
        })))
        geneMean = mean(apply(matr2,1,sum)) 
        if(geneMean<gene_base_mean_thresh) return(NULL)
        baseMean=apply(matr2,2,mean) ## can use matr2 or matr1
         if(length(which(baseMean>0))<min_num_transcripts) return(NULL)
         inds_t =which(baseMean>transcript_base_mean_thresh)
         if(length(inds_t)==0) return (NULL)
         totsum1 = apply(matr1,1,sum)
        totsum2 = apply(matr2,1,sum)
        if(sum(totsum2)<min_total_sum) return(NULL)
        matr1 = cbind(totsum1, matr1[,inds_t, drop=F])
        matr2 = cbind(totsum2, matr2[,inds_t, drop=F])
        dimnames(matr1)[[2]][1]=attr(matr,"gene")
        dimnames(matr2)[[2]][1]=attr(matr,"gene")
        
        ##WORKING HERE
        if(length(psb1$def$adjust_fixed)>0){
          tobind = psb1$df_all[,match(psb1$def$adjust_fixed,names(psb1$df_all)),drop=F]
          for(k in 1:ncol(tobind)) tobind[,k] = as.numeric(tobind[,k])
         matr1 = cbind(tobind, matr1)
         tobind = psb1$df_summ[,match(psb1$def$adjust_fixed,names(psb1$df_summ)),drop=F]
         for(k in 1:ncol(tobind)) tobind[,k] = as.numeric(tobind[,k])
         matr2 = cbind(tobind, matr2)
        }
#        means_all = c(baseMean[inds_t], geneMean)
#        c(gene,dimnames(matr1)[[2]]), attr(psb1$def,"nme"),Condition,psb1$cell_type, 
        
        pv_11=.assocTestAll(psb1, matr1, matr2, adj = 1:(1+length(psb1$def$adjust_fixed)), var_thresh = var_thresh)
        attr(pv_11,"nme")=attr(psb1$def,"nme")
      pv_11
      })
      ##SOMETHING WRONG WITH MERGING and numeric columns
      
     results_all = .merge1(results1, all_num=T,addName="Cell",addRow="ID",
                           attribute = list(Name="nme"))
     if(is.null(results_all)) return (NULL)
     if(nrow(results_all)==0) return (NULL)
     
     results_all[match(header1, names(results_all))]
  })
  
  for(j in 1:length(outf)){
         write.table(pv_res_all[[j]],outf[[j]],row.names=F,quote=F, sep="\t", append=T, col.names=F)
  }

}
                          
closeAllConnections()
}
#.runType("Cell_State" )



