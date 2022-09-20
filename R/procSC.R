if(FALSE){

  setwd( "/data/scratch/projects/punim1068/analysis_genome/resdir8")
  setwd("/scratch/punim1068/flames_analysis")
  #setwd("C:/Users/LCOIN/Downloads/scp")
#setwd("/scratch/punim1068/Sepsis_single_cell")
  setwd("C:/Users/LCOIN/data/flames_multigo")
  assign("last.warning", NULL, envir = baseenv())
  rm(list = ls())
  closeAllConnections()
  
#load(".RData")
}
  
library(lmerTest)
library(jsonlite)

params_file="params.json"
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0) params_file = args[1]


print(params_file)
params = read_json(params_file); options(params)  
outdir = params$outdir
if(length(args)>1) {
 outdir = args[2]
}else{
 outdir=paste0(params$outdir,format(Sys.time(), "%a_%b_%d_%H_%M_%S_%Y"))
}
params$outdir=outdir
working_dir = getOption("working_dir","./")
setwd(working_dir)
#try(source("/data/gpfs/projects/punim1466/sw/multigo/R/procSC_lib.R"))
try(source(paste(getOption("procSC_home","/home/lcoin/github/procSC"),"R/procSC_lib.R",sep="/")))

split = params$split[[1]]
collapse=params$collapse[[1]]


samples=read.csv(params$meta_file[[1]],sep="\t",header=T)
##NOW INTEGRATE ANY CLINICAL DATA WITH SAMPLE TABLE
samples=.addClinicalData(samples, params$clinical_data)
##FOLLOWING IF WE DONT HAVE TOTAL READS COUNTS AS COLUMN IN meta_file
if(nchar(getOption("sumsInMeta",""))==0 || is.null(samples[getOption("sumsInMeta")])){
  samples=.addCountData(samples)
}
samples=.addSplitColumns(samples, params$cohortSplitName,params$cohortSplit)
samples=.addMergeColumns(samples, getOption("mergeColumns"), params$cohortSplit[[1]])



#def_all=list( splitBy="all", case=list(Cohort="all"), control = list(Cohort="none"))
#subindices_samples=.getSubIndices(def_all, samples,min_cells_cases=0 , min_cells_controls=0)
#c(list(all=def_all),

defs_=.readDefs(params$defs)

#samp2 = lapply(defs_, .splitFrameMultiple,samples1)
##NOW READ SC data

connection = gzfile( params$input_file[[1]],"rb")
header=strsplit(readLines(connection,1),params$split[[1]])[[1]]
samples1 = lapply(defs_, function(def){
  samp1 = .getDataFrame(def, samples, header)
  tab = table(samp1$casecontrol)
  print(min(tab))
  samp2 = .splitFrameMultiple(def, samp1)
  
  psB=lapply(samp2,.pseudoBulk, params$pseudobulk,def,  sumTags = params$sumsInMeta)
  print(lapply(psB, function(x)table(x$df_all$casecontrol)))
  #attr(psB, "def")=def
  psB
})

#psb1=samples1$VIC_infected_vs_uninfected_child$Uninfected
dir.create(outdir)

outf=vector("list", length(samples1))

for(j in 1:length(samples1)){
  outf[[j]] = gzfile(paste(outdir, "/",names(samples1)[[j]],".tsv.gz",sep=""),open="w")
}

nextLine=strsplit(readLines(connection,1),params$split[[1]])[[1]]


while(TRUE){
  matr=.readMatrix(connection,nextLine)
  nextLine = attr(matr,"next")
   
  pv_res_all =  lapply(samples1,function(pSB){
     params_ = attr(psB,"params")
      lapply(pSB, function(psb1){
      .assocTestAll(psb1, matr, ll0)
      })
  })
  for(j in 1:length(outf)){
    lapply(pv_res_all[[j]], function(pv_res){
      write.table(pv_res,outf[[j]],row.names=F,quote=F, sep="\t")
    })
  }
 
}
                          
closeAllConnections()
#.runType("Cell_State" )



