if(FALSE){

  setwd( "/data/scratch/projects/punim1068/analysis_genome/resdir8")
  setwd("/scratch/punim1068/flames_analysis")
  #setwd("C:/Users/LCOIN/Downloads/scp")
#setwd("/scratch/punim1068/Sepsis_single_cell")
  assign("last.warning", NULL, envir = baseenv())
  rm(list = ls())
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

#try(source("/data/gpfs/projects/punim1466/sw/multigo/R/procSC_lib.R"))
try(source("/home/lcoin/github/procSC/R/procSC_lib.R"))

split = params$split[[1]]
collapse=params$collapse[[1]]

closeAllConnections()

samples=read.csv(params$meta_file[[1]],sep="\t",header=T)

##NOW INTEGRATE ANY CLINICAL DATA WITH SAMPLE TABLE
if(!is.null(params$clinical_data)){
  clinical_data = read.table(params$clinical_data[[1]], sep="\t", header=T)
  nv = rep(0, ncol(clinical_data))
  for(j in 1:ncol(clinical_data)){
    nv[j] =length(unique(clinical_data[,j]))
  }
  clinical_data = clinical_data[,nv>1]  
  mi = match( samples$Patient,clinical_data$Patient )
  nonNA = which(!is.na(mi))
  samples = samples[nonNA,]
  mi = mi[nonNA]
  clinical_data1 = array(dim =c(nrow(samples),ncol(clinical_data)), dimnames = list(dimnames(samples)[[1]], names(clinical_data)))
  
  mi = match( samples$Patient,clinical_data$Patient )
  for(j in 1:nrow(clinical_data)){
    j1 = which(mi==j)
   # j1 = j1[!is.na(j1)]
    if(length(j1)>0){
      for(i in 1:ncol(clinical_data)) clinical_data1[j1,i] = clinical_data[j,i]
    }
  }
  samples = cbind(samples, clinical_data1[,-1])
}

cohortSplit = params$cohortSplit[[1]]
if(!is.null(cohortSplit)){
if(nchar(cohortSplit)>0){
  nmes1 = strsplit(params$cohortSplitName[[1]],":")[[1]]
 cohort1= lapply(as.character(samples$Cohort), function(x){
   spl=strsplit(x,cohortSplit)[[1]]
   diff1 = length(nmes1) - length(spl)
   if(diff1>0) return(c(rep("NA", diff1),spl))
   else return(spl)
 })
 df1 = t(data.frame(cohort1))
 dimnames(df1) = list(NULL,nmes1)
 samples =cbind(samples, df1)
}
}

if(!is.null(getOption("mergeColumns"))){
  mergeColumns = getOption("mergeColumns")
  new_samples = data.frame(matrix(nrow = nrow(samples), ncol = length(mergeColumns)))
  names(new_samples) = names(mergeColumns)
  for(k in 1:length(mergeColumns)){
    new_samples[,k] = apply(samples[,match(mergeColumns[[k]], names(samples) )],1,paste,collapse=cohortSplit)
  }
  samples = cbind(samples, new_samples)
}



connection = gzfile( params$input_file[[1]],"rb")
header=strsplit(readLines(connection,1),params$split[[1]])[[1]]
mi = match(header, samples$NAME)
##which barcodes have been removed
nonNA_mi = which(!is.na(mi))
if(nchar(getOption("sumsInMeta",""))>0){
 sums= samples[[getOption("sumsInMeta")]][mi]
 names(sums) = header
}else{
  sums = .getSums()
  names(sums) = header
  ##following makes a new meta file with the counts
  outf1 = paste0(params$meta_file,".1.txt")
  if(!file.exists(outf1)){
    mi1 = match( samples$NAME,header)
    Count = sums[mi1]
   # indsNA1 = which(is.na(names(Count)))
    names(Count) = dimnames(samples)[[1]]
    samples1 = cbind(samples, Count)
    write.table(samples1, file=outf1,sep="\t",quote=F, row.names=F)
    samples = samples1
  }
}
if(length(nonNA_mi)==0) stop(" cannot match")
print(paste("total matching cell " ,length(nonNA_mi)))

all = rep("all", length(mi))
##samples 1 re ordered to match header
samples1 = cbind(samples[mi,],all)

##CHECKING THINGS MATCHED UP PROPERLY
chk=cbind(header[nonNA_mi],as.character(samples1$NAME[nonNA_mi]))
if(length(which(apply(chk,1,function(x)x[1]!=x[2])))>0){
  stop("error with matching")
}


for(j in 2:ncol(samples1)) samples1[,j] = factor(samples1[,j])

def=list( splitType="all", type="Cohort",case= unique(samples1$Cohort),control = c())
subindices_samples=.getSubIndices(def, samples1,min_cells_cases=0 , min_cells_controls=0)
                                 

all_sum=.processLine(sums, rep(1,length(sums)), subindices_samples$all,gene_name = "all", debug=F, canSkip=F)
sums_all = rep(NA, length(sums))
for(i in 1:length(all_sum$sums_cases)){
  indsi = subindices_samples$all$cases[subindices_samples$all$inds_cases[[i]]]
  sums_all[indsi] = all_sum$sums_cases[i]
}
names(sums_all) = names(sums)

meta=.getMeta(samples1)


#defs_all2= rev(defs_all2)
#options(params)
if(is.null(names(params$defs))) names(params$defs)=1:length(params$defs)

for(i in 1:length(params$defs)){
  print(i)
  defs_ = .readDefs(params$defs[[i]])
  for(j in 1:length(defs_)){
    def_nme =names(defs_)[[j]]
   .runType(params, defs=defs_[[j]], def_nme=def_nme, outdir=outdir)
  }
}

print(" done types")
closeAllConnections()
#.runType("Cell_State" )



