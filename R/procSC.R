if(FALSE){
  assign("last.warning", NULL, envir = baseenv())
  rm(list = ls())
  closeAllConnections()
  setwd("~/github/procSC/example/flames_jessie")
}
{  
library(lmerTest)
library(jsonlite)

params_file="params.json"
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) params_file = args[1]


print(params_file)
params = read_json(params_file); options(params)  
if(!is.null(getOption("include"))){
  params1 = read_json(getOption("include")); options(params1); options(params);
}
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
defs_=.readDefs(getOption("defs","defs.json"))

whitelist = NULL
if(!is.null(getOption("transcript_whitelist"))){
  whitelist = read.csv(getOption("transcript_whitelist"),head=F)[,1]
}

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

print(head(samples))

connection = gzfile( getOption("input_file"),"rb")
header=strsplit(readLines(connection,1),getOption("split",",")[[1]])[[1]]

samples1 = lapply(defs_, function(def){
  mi = match( names(def$case),names(samples))
  names(mi)=names(def$case)
  if(length(which(is.na(mi)))>0) stop("paste problem with defs ")
  samp1 = .getDataFrame(def, samples, header)
  if(nrow(samp1)==0) stop(" somthing wrong with definition, no barcodes selected")
  tab = table(samp1$casecontrol)
  print(tab)
  samp2 = .splitFrameMultiple(def, samp1)
  psB=lapply(samp2,.pseudoBulk, def)
           
  print(lapply(psB, function(x)table(x$df_all$casecontrol)))
  #attr(psB, "def")=def
  psB
})

#psb1=samples1$VIC_infected_vs_uninfected_child$all
dir.create(outdir)

outf=vector("list", length(samples1))
header =  c("ID","Name","Condition","Cell_type","pv_cell","beta_cell","pv_pb","beta_pb")
for(j in 1:length(samples1)){
  outf[[j]] = gzfile(paste(outdir, "/",names(samples1)[[j]],".tsv.gz",sep=""),open="w")
  writeLines(paste(header,collapse="\t" ), outf[[j]])
}

nextLine=strsplit(readLines(connection,1),getOption("split")[[1]])[[1]]

cnt=0
max_count=getOption("max_cnt",1e9)
min_num_reads = getOption("min_num_reads")

transcript_column = getOption("transcriptColumn",NULL)
while(cnt<max_count){
  print(paste(cnt, nextLine[[1]]))
  
  matr=.readMatrix(connection,nextLine)
  nextLine = attr(matr,"next")
  
  to_include =rep(T, length(matr))
  if(!is.null(transcript_column) && !is.null(whitelist)){
    transcript_ids = unlist(lapply(names(matr), function(m1){
      strsplit(m1,"\\|")[[1]][transcript_column]
    }))
    to_include = (!is.na(match(transcript_ids, whitelist)))
  }
   
  pv_res_all =  lapply(samples1,function(psB){
     params_ = attr(psB,"params")
      lapply(psB, function(psb1){
      .assocTestAll(psb1, matr, to_include=to_include, min_num_reads = min_num_reads)
      })
  })
  for(j in 1:length(outf)){
    lapply(pv_res_all[[j]], function(pv_res){
         write.table(pv_res,outf[[j]],row.names=F,quote=F, sep="\t", append=T, col.names=F)
    })
  }
 cnt=cnt+1
}
                          
closeAllConnections()
}
#.runType("Cell_State" )



