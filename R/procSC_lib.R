
#
.sumLines<-function(lines,split=getOption("split","\t")){
  splitlines = lapply(lines, function(line)strsplit(line,split)[[1]] )
  sum = rep(0, length(splitlines[[1]]))
  for(i in 1:length(splitlines)){
    inds = which(splitlines[[i]]!="0")[-1]
    sum[inds] = sum[inds]+ as.numeric(splitlines[[i]][inds])
  }
  sum
}

readLines2<-function(connection, n,nxt, normaliseByGene = getOption("normaliseByGene",FALSE),
                     geneColumn = getOption("geneColumn",1), split = getOption("split","\t")){
  line=readLines(connection,n)
  if(length(line)==0) {
    if(is.null(nxt$line)) return(NULL)
    lines = list(nxt$line)
    result = list(lines = lines, line=NULL, gene=NULL)
    return(result)
  }
  l1=substr(line,1,gregexpr(pattern=split ,line)[[1]][1]-1)#,"|")
  gene=strsplit(l1,"\\|")[[1]][geneColumn]
  if(is.null(nxt)){
    result = list(line=line, gene=gene)
  }else{
    gene0 =nxt$gene
    lines = list(nxt$line)
    
    while(gene==gene0){
       lines = c(lines,line)
       line=readLines(connection,n)
       l1=substr(line,1,gregexpr(pattern=split ,line)[[1]][1]-1)#,"|")
       gene=strsplit(l1,"\\|")[[1]][geneColumn]
    }
  
    result = list(lines = lines, line=line, gene=gene)
  }
  return(result)
}

readLines1<-function(connection, n,normaliseByGene = getOption("normaliseByGene",FALSE),
                     geneColumn = getOption("geneColumn",0), split = getOption("split","\t")){
  line =readLines(connection,n)
  if(length(line)>0 && !is.null(getOption("subFirst",NULL))){
    if(getOption("subFirst", FALSE)[[1]]) line = sub(getOption("split",",")[[1]],"|",line)
  }
 
  line
}
##FUNCTIONS
.restrict<-function(samples1, restrict, inds_r){
 if(length(restrict)==0) return(inds_r)
 for(k in 1:length(restrict)){
     inds_r = inds_r & samples1[[names(restrict)[k]]] %in% restrict[[k]] 
    }
inds_r
}


#.getInds<-function(def,samples1,patient, header_out, outdir, collapse,meta, cell_types){

#  .getSubIndices(def, samples1, patient, header_out, outdir1, "", collapse, def$restrict)
#}

.simplify<-function(pc,sep="_", simplify=F){
  if(!simplify) return(factor(pc))
  m1=t(data.frame(lapply(pc, function(x) strsplit(x,sep)[[1]])))
  l1=apply(m1,2, function(v) length(unique(v)))
  factor(apply(m1[,l1>1,drop=F], 1,function(v) paste(v, collapse=sep)))
}
#subindices_type = lapply(defs, .getSubIndices, samples1, header_out, outdir=outdir_nme, outf_cell = outf_cell,
#                         collapse=collapse, meta=samples_p2, cell_types=cell_types)
#subindices_samples=.getSubIndices(def, samples1,min_cells_cases=0 , min_cells_controls=0)


.getSubIndices<-function(def, samples1,  
                        header_out = "", outdir=NULL, outf_cell = NULL,
                         collapse="\t", meta=NULL, cell_types=NULL,
                        min_cells_cases = getOption("min_cells_cases",1), 
                        min_cells_controls = getOption("min_cells_controls",1)
){
  pv_index = getOption("pvs_index",11);pv_thresh = getOption("pvs_thresh",0.01);out_nme=attr(def,"nme")
  patients = samples1$Patient
  lmTest = getOption("lmTest","combined")
  if(is.null(def$splitType)) stop("error - did not define split type")
  if(!is.null(outdir) && !is.null(meta)){
    dir.create(outdir)
   outdir1 = paste(outdir, out_nme,sep="/")
   dir.create(outdir1)
   write.csv(meta, paste(outdir1,"Metadata_scp.txt",sep="/"))  
   write.table(paste0(paste0("_",cell_types),"."), paste(outdir1,"cell_types.txt",sep="/"), quote=F, row.names=F, col.names=F,sep=",")
  }
  type =   def$splitType; col_nme = def$type;sepsis =  def$case; control=def$control;restrict=def$restrict 
  ##THIS ALLOWS FOR COMPOUND TYPES  
  print(paste("TYPE", type))
  if(length(which(names(samples1)%in% col_nme))==0){
    stop("problem")
  }
  if(length(which(names(samples1) %in% type))==0){
     samples1 = cbind(samples1,
                      apply( samples1[,which(names(samples1) %in% strsplit(type,",")[[1]]),drop=F],1,paste,collapse="_"))
     names(samples1)[ncol(samples1)]=type

  }

   print(paste("col_nme", col_nme))
   cell_types = unique(samples1[[type]])
  cohorts = grep("group",levels(samples1[[col_nme]]),  inv=T,v=T)
  if(is.null(control)){
    control = cohorts[!(cohorts %in% sepsis)]
  }
  if(is.null(sepsis)){
    sepsis = cohorts[!(cohorts %in% control)]
  }
  adjust=unlist(def$adjust)
  
  subindices=lapply(cell_types,function(cell){
    cases = which(samples1[[type]] %in% cell & samples1[[col_nme]] %in% sepsis )
    controls =  which(samples1[[type]] %in% cell & samples1[[col_nme]] %in% control)
    case_ids = unique(patients[cases])
    control_ids = unique(patients[controls])
    df_cases = data.frame(matrix(NA,nrow=length(cases), ncol=4))
    df_controls =  data.frame(matrix(NA,nrow=length(controls), ncol=4))
    pcomb_cases = .simplify(patients[cases])
    pcomb_controls =.simplify(patients[controls])
    comb_inds_cases=lapply(levels(pcomb_cases), function(y){
       which(pcomb_cases %in% y)
    })
    comb_inds_controls=lapply(levels(pcomb_controls), function(y){
      which(pcomb_controls %in% y)
    })
    names(comb_inds_cases) =levels(pcomb_cases)
    names(comb_inds_controls) = levels(pcomb_controls)
    case_ids = names(comb_inds_cases)
    control_ids = names(comb_inds_controls)
    if(length(comb_inds_cases)>0){
      for(jk in 1:length(comb_inds_cases)){
        df_cases[comb_inds_cases[[jk]],2] = names(comb_inds_cases)[jk]
        df_cases[comb_inds_cases[[jk]],3] =cases[comb_inds_cases[[jk]]]
      }
    }
    if(length(comb_inds_controls)>0){
      for(jk in 1:length(comb_inds_controls)){
        df_controls[comb_inds_controls[[jk]],2] = names(comb_inds_controls)[jk]
        df_controls[comb_inds_controls[[jk]],3] =controls[comb_inds_controls[[jk]]]
      }
    }
    df_cases[,1]=rep(1, nrow(df_cases)); df_controls[,1] = rep(0, nrow(df_controls))
    df = rbind(df_cases, df_controls)
   
    names(df) = c("condition","replicate","index","value")
    if(!is.null(adjust)) {
      df_to_add = samples1[c(cases, controls), match( adjust, names(samples1)),drop=F]
      names(df_to_add) =adjust 
      df  = cbind(df,df_to_add)
    }
    samps_cases=unlist(lapply(comb_inds_cases, length))
    samps_controls=unlist(lapply(comb_inds_controls, length))
    tsum_x = sum(samps_cases,na.rm=T) 
    tsum_y = sum(samps_controls, na.rm=T)
    df1_cases = data.frame(matrix(NA,nrow=length(samps_cases), ncol=3))
    df1_controls =  data.frame(matrix(NA,nrow=length(samps_controls), ncol=3))
    names(df1_cases) = c("condition","replicate","value")
    names(df1_controls) = c("condition","replicate","value")
    df1_cases$condition=rep(1, nrow(df1_cases));  df1_controls$condition = rep(0, nrow(df1_controls))
    df1_cases$replicate = names(samps_cases)
    df1_controls$replicate = names(samps_controls)
    if(length(samps_cases)>0) df1_cases$value = samps_cases
    if(length(samps_controls)>0) df1_controls$value = samps_controls
    df1 = rbind(df1_cases, df1_controls)
  
    df[,1] = factor(df[,1]); df[,2] = factor(df[,2]) ;
    df1[,1] = factor(df1[,1]); df1[,2] = factor(df1[,2]) ;
    
    if(!is.null(adjust)){
      df1_to_add = meta[match(df1$replicate, meta$Patient), match( adjust, names(meta)),drop=F]
     names(df1_to_add) = adjust
     df1 = cbind(df1, df1_to_add)
    }
      pv=.wilcoxon(x = samps_cases,y=  samps_controls,
                   paste("cell_count",cell,sep="_"),tsum_x, tsum_y)
      pvi = as.numeric(pv[pv_index])
      if(lmTest=="combined" && nrow(df1)>2){
        if(!is.na(pvi) && pvi<pv_thresh){
        pv1=lmtestAnalysis(df1,paste("cell_count",cell,sep="_"),tsum_x, tsum_y)
        pv[pv_index] = pv1[pv_index]  
      }
     }
    outf=NULL
    cell_type=cell
    if(!is.null(outf_cell))writeLines(paste(c(pv,as.character(cell),out_nme), collapse=collapse),outf_cell)
    if(length(cases)<min_cells_cases || length(controls) < min_cells_controls){
       print(paste("excluding ",out_nme,cell,length(cases), length(controls), min_cells_cases, min_cells_controls))
       return (NULL)
    }
    if(!is.null(outdir)){
      fi = paste(paste(outdir1, paste(out_nme,cell,sep="_"), sep="/"),"csv",sep=".")
      print(paste("opening ",fi))
        outf = file(fi, open="w") 
      writeLines(header_out, outf,sep=collapse)
      writeLines(paste(c(names(comb_inds_cases), names(comb_inds_controls)), collapse=collapse), outf)
      writeLines(paste(c(pv,as.character(cell)), collapse=collapse),outf, sep=collapse)
      writeLines(paste(c(samps_cases, samps_controls), collapse=collapse), outf)
    }
    list(cases=cases,controls=controls, inds_cases = comb_inds_cases,inds_controls = comb_inds_controls,
         nme=out_nme, control_ids = control_ids, case_ids = case_ids, outf=outf, cell_type = cell_type, df = df, df1 = df1)
          })
  names(subindices)=cell_types
  
  subindices[!unlist(lapply(subindices, function(x) is.null(x)))]
}


lmtestAnalysis<-function(df1,gene="", tsum_x=10, tsum_y=10,  lmer=T){
 # x$replicate = factor(x$replicate)
  inds_cases = which(df1$condition==1 & !is.na(df1$value))
  inds_controls = which(df1$condition==0& !is.na(df1$value))
  lx=length(inds_cases)
  ly = length(inds_controls)
  sx = sum(df1$value[inds_cases])
  sy = sum(df1$value[inds_controls])
  mx = sx/lx
  my = sy/ly
  baseMean = (sx+sy)/(lx+ly)
  pvalue=NA
  beta=NA
  if(lx>=2 && ly>=2){
    if(lmer){
      o = try(lmerTest::lmer(value ~ condition + (1 | replicate), data=df1))
    }else{
      o = try(lmerTest::lm(df1$tail_length ~ df1$condition))
    }
    
    if(!inherits(o,"try-error")) {
      coeff = summary(o)$coeff
      beta=coeff[2,1]
      pvalue = coeff[2,5]
    }
  } 
#  c(gene, baseMean,  mx, my,beta ,lx, ly,tsum_x, tsum_y, tsum_x+tsum_y, pvalue)
  c(gene,sprintf("%5.3g",c(baseMean,  mx, my,beta ,lx, ly,tsum_x, tsum_y, tsum_x+tsum_y, pvalue)))
}

##min 5 reads
.wilcoxon<-function(x,y,gene="", tsum_x=10, tsum_y = 10){
  nonNAx = which(!is.na(x))
  nonNAy = which(!is.na(y))
  lx= length(nonNAx)
  ly = length(nonNAy)
  sx = sum(x[nonNAx])
  sy = sum(y[nonNAy])
  mx = if(lx==0) NA else median(x[nonNAx])
  my = if(ly==0) NA else  median(y[nonNAy])
  baseMean = (sx+sy)/(lx+ly)
 # if(is.infinite(mx)) stop(paste("inf mx", sx, lx))
 # if(is.infinite(my)) stop(paste("inf y", sy, ly))
  if(lx>=2 && ly>=2){
   pvalue= try(wilcox.test(x[nonNAx],y[nonNAy])$p.value)
   if(inherits(pvalue,"try-error")) pvalue=NA
   } else pvalue= NA
  c(gene,sprintf("%5.3g",c(baseMean,  mx, my,log2(mx/my) ,lx, ly,tsum_x, tsum_y, tsum_x+tsum_y, pvalue)))
}


##nonNA_mi is the subindices to consider
#  .processLine(sums,sums, subinds$inds, subinds$out_files,splitline[1], normalise=T)
#.processLine(sums,sums, subinds,"totalSum", normalise=T)
#


.processLine<-function(splitline,sums1, subinds, gene_name="",sample_sums=NULL, collapse="\t", 
                       debug=F,lmTest=getOption("lmTest","wilcoxon"), canSkip=T, outf_cell = NULL,
                       pv_index = getOption("pvs_index",11),
                       pv_thresh = getOption("pvs_thresh",0.01),
                       min_num_reads = getOption("min_num_reads",10)){
  mult=1
  x=(splitline[subinds$cases])
  y=(splitline[subinds$controls])
      tsum_x = sum(x, na.rm=T) 
      tsum_y = sum(y, na.rm=T)
      if(canSkip && tsum_x+tsum_y < min_num_reads) return (NULL)
    if(lmTest=="lmTest"){
      df = subinds$df
      df$value = (splitline[df$index])/sums1[df$index]
      pvs=lmtestAnalysis(df,gene_name,tsum_x, tsum_y)
      
      df1=subinds$df1
      df1$value= c( unlist(lapply(subinds$inds_cases, function (v){
          ind1 = subinds$cases[v]
          sum(splitline[ind1])/sum(sums1[ind1])
        })),
        unlist(lapply(subinds$inds_controls, function (v){
          ind1 = subinds$controls[v]
          sum(splitline[ind1])/sum(sums1[ind1])
        })))
      pvs1 = lmtestAnalysis(df1, gene_name, tsum_x, tsum_y)
    }else{
      x_sum =sums1[subinds$cases]
      y_sum=sums1[subinds$controls]
      x_=mult*(x/x_sum)
      y_=mult*(y/y_sum)
      
      pvs =  .wilcoxon(x_,y_, gene_name,tsum_x, tsum_y)
      pvi = as.numeric(pvs[pv_index])
      
      x1=unlist(lapply(subinds$inds_cases, function(v) sum(x[v])))
      y1=unlist(lapply(subinds$inds_controls, function(v) sum(y[v])))
      x1_sum=unlist(lapply(subinds$inds_cases, function(v) sum(x_sum[v])))
      y1_sum=unlist(lapply(subinds$inds_controls, function(v) sum(y_sum[v])))
      x1_=mult*(x1/x1_sum)
      y1_=mult*(y1/y1_sum)
      pvs1 = .wilcoxon(x1_,y1_, gene_name,tsum_x, tsum_y)
    #  pvi1 = as.numeric(pvs1[pv_index])
#    print(paste("here", pvi, lmTest, as.numeric(pvi)))
    if(lmTest=="combined"){
      if(!is.na(pvi) && pvi<pv_thresh){
        df = subinds$df
        df$value = (splitline[df$index])/sums1[df$index]
        pvs_1=lmtestAnalysis(df,gene_name,tsum_x, tsum_y)
        pvs[pv_index] = pvs_1[pv_index]
        
        df1=subinds$df1
        df1$value= c( unlist(lapply(subinds$inds_cases, function (v){
          ind1 = subinds$cases[v]
          sum(splitline[ind1])/sum(sums1[ind1])
        })),
        unlist(lapply(subinds$inds_controls, function (v){
          ind1 = subinds$controls[v]
          sum(splitline[ind1])/sum(sums1[ind1])
        })))
        pvs1_1=lmtestAnalysis(df1,gene_name,tsum_x, tsum_y)
        pvs1[pv_index] =pvs1_1[pv_index]
       # print(paste("lm test",pvi,"-->",pvs[pv_index]))
        
      }
    }
    }
	if(debug){
		names(pvs)= c("GeneName","baseMean", "median_x","median_y","log2FoldChange", "length_x","length_y","tsum_x","tsum_y","tsum", "pvalue")
		 infx = which(x_>1)
		 infy = which(y_>1)
		    if(length(infx)>0){
			print("is infinite")
			print(cbind(x,sums1[subinds$cases] )[infx,])
			stop(paste(gene_name, "END"))
			}
if(length(infy)>0){
			print("is infinite")

			print(cbind(y,sums1[subinds$controls])[infy,])
			stop(paste(gene_name, "END"))
			}	
	
	}
      pvi = as.numeric(pvs[pv_index])
      if(is.na(pvi)  && canSkip){
       return( list(pvs=pvs))
      }
  sums_cases = unlist(lapply(subinds$inds_cases, function(x1){
    sum(x[x1], na.rm=T)
  }))
  sums_controls = unlist(lapply(subinds$inds_controls, function(x1){
    sum(y[x1], na.rm=T)
  }))
  if(!is.null(outf_cell)){
    writeLines(paste(pvs, collapse=collapse), outf_cell,sep=collapse)
    writeLines(paste(pvs1, collapse=collapse), outf_cell,sep=collapse)
    writeLines(as.character(subinds$cell_type), outf_cell,sep="\t")
    writeLines(as.character(subinds$nme), outf_cell,sep="\n")
  }
  if(!is.null(subinds$outf)){
    #print(paste(sum(sums_cases)== sum(x), sum(sums_controls)== sum(y)))
    
    writeLines(paste(pvs, collapse=collapse), subinds$outf,sep=collapse)
    writeLines(paste(pvs1, collapse=collapse), subinds$outf,sep=collapse)
    
      writeLines(as.character(subinds$cell_type), subinds$outf,sep=collapse)
     
    writeLines(paste(sums_cases, collapse=collapse), subinds$outf,sep=collapse)
    writeLines(paste(sums_controls, collapse=collapse), subinds$outf,sep="\n")
    if(!is.null(sample_sums)){
      pvs_1 =rep("NA", length(pvs))
      pvs_1[1] = "sampleSums"
      writeLines(paste(pvs_1, collapse=collapse), subinds$outf,sep=collapse)
      writeLines(paste(pvs_1, collapse=collapse), subinds$outf,sep=collapse)
      writeLines(as.character(subinds$cell_type), subinds$outf,sep=collapse)
      writeLines(paste( sample_sums[match(names(sums_cases), names(sample_sums))], collapse=collapse), subinds$outf,sep=collapse)
      writeLines(paste( sample_sums[match(names(sums_controls), names(sample_sums))], collapse=collapse), subinds$outf,sep="\n")
    }
  }

  if(debug){
	res1 = list(cases = x, controls = y,norm_cases = sums1[subinds$cases], norm_controls = sums1[subinds$controls], pvs=pvs)
	return(res1)
  }
 list(pvs=pvs,pvs1=pvs1,sums_cases=sums_cases, sums_controls=sums_controls)
}







.getMeta<-function(samples1){
  samples_t = table(samples1$Patient)                                       
  samples_p1 = samples1[!duplicated(samples1$Patient),,drop=F]
  mi1 = match(samples_p1$Patient, names(samples_t))
  cell_counts = samples_t[mi1]
  samples_p2 = cbind(samples_p1, cell_counts)
  toremove = names(samples_p2) %in% c("NAME","Cell_Type", "Cell_State")
  print(toremove)
  samples_p2 = samples_p2[,!toremove,drop=F]
  samples_p2
}
##OUTPUT FOR TYPE
.runType<-function(params, defs, def_nme, outdir=getOption("outdir")){
  min_sum=getOption("min_total_sum",100)
  geneColumn=getOption("geneColumn"); input_file=getOption("input_file"); 
 #sum_file=getOption("sum_file",paste0("sum_", input_file))
 normalise= getOption("normalise");split=getOption("split"); collapse=getOption("collapse"); lmTest = getOption("lmTest")
  closeAllConnections()
  h_out1 = c("GeneName","baseMean", "mean_x","mean_y","log2FoldChange", "length_x","length_y","tsum_x","tsum_y", "tsum",  "pvalue")
  if(lmTest=="lmTest"){
    h_out1[5] = "coeff"
  }
   h_out = c(h_out1, paste(h_out1,"pb",sep="."), "cell_type")
  
  
  header_out= paste(h_out,collapse=collapse)
  pv_index = which(h_out=="pvalue")
  dir.create(outdir)
  outdir_nme =paste(outdir,def_nme,sep="/") 
  dir.create(outdir_nme)
  write_json(params, paste(outdir,"params_out.json",sep="/"), pretty=T, auto_unbox=T)
  write_json(defs,paste(outdir_nme,"defs_out.json",sep="/"), pretty=T, auto_unbox=T)
  
  #samples_p2$Patient = apply(cbind(samples_p2$Patient, samples_p2$Cohort),1,paste,collapse="_")
  cell_types = as.matrix(grep("group",unique(samples1$Cell_Type),inv=T,v=T), ncol=1)
  cell_types = cell_types[!is.na(cell_types)]
  outf_cell = file(paste0(outdir_nme,"/cell_sig.txt"), open="w") 
  writeLines(paste(header_out,"nme",sep=collapse), outf_cell)
  
  subindices_type = lapply(defs, .getSubIndices, samples1, header_out, outdir=outdir_nme, outf_cell = outf_cell,
                           collapse=collapse, meta=meta, cell_types=cell_types)
  names(subindices_type) = names(defs)
  subindices1 = unlist(subindices_type, recursive=F)
if(length(subindices1)==0) stop(" no comparisons have sufficient cell counts")
  pvs = lapply(subindices1, function(subinds){
    .processLine(sums,sums_all, subinds,"totalSum",sample_sums = all_sum$sums_cases, outf_cell = outf_cell,
                 collapse=collapse, canSkip=F)
  })
  close(outf_cell)
  outf_cell = NULL
  connection = gzfile( input_file,"rb")
  connection1 = NULL
  if(getOption("normaliseByGene",FALSE)){
    #connection1 reads the gene counts in same order or transcript counts
   connection1 =   gzfile( sum_file,"rb")
   header1=strsplit(readLines(connection1,1),split)[[1]]
  # print(substr(h1,1,30))
  }
  header=strsplit(readLines(connection,1),split)[[1]]
 
  if(!normalise){
    sums_1 = rep(1, length(sums))
  }else{
    sums_1 = sums
  }
  line = readLines(connection,1) #
  gene_counts = NULL
  gene1= "" 
  
  k=0
  k1=0
  max_lines = getOption("max_lines",3e6) 
  
  tryCatch(
    {
  while(!is.null(line) && length(line)>0 && k<max_lines){
    splitline = strsplit(line,split[[1]])[[1]]
    line = readLines(connection,1) # next line
    k = k+1

    inds2_=which(splitline!="0")
    splitline1 = rep(0,length(splitline))
    splitline1[inds2_] = as.numeric(splitline[inds2_])
    sum_inds2 = sum( splitline1[inds2_], na.rm=T)
    if(sum_inds2<min_sum){
      next;
    }
    if(getOption("normaliseByGene",FALSE)){
       ##need something here
       gene = strsplit(splitline[[1]],"\\|")[[1]][getOption("geneColumn",0)]
       if(gene!=gene1){
         line1 = readLines(connection1,1) 
         gene_counts = strsplit(line1,split)[[1]]
         gene1 = gene_counts[1]
         k1=k1+1
         jj=1
         while(gene1!=gene){
           line1 = readLines(connection1,1) 
           gene_counts = strsplit(line1,split)[[1]]
           gene1 = gene_counts[1]
           jj = jj+1
         #  print(jj)
         }
        # if(gene!=gene1) stop(" the gene counts and transcript counts are not in same order")
         print(paste(gene,k1))
         inds2_gene=which(gene_counts!="0")
         sums_1 = rep(0, length(gene_counts))
         sums_1[inds2_gene] = as.numeric(gene_counts[inds2_gene])
#         maxdiff = max(sums_1[inds2_gene] - as.numeric(splitline[inds2_gene]),na.rm=T)
       }
    }
    tryCatch(
      {
    	if(params$debug[[1]]){
    		merge = cbind(as.numeric(splitline[-1]), sums_1[-1])
    		inds1 = which(apply(merge, 1, function(v) v[1]>v[2]))
    		if(length(inds1)>0) stop("error")
    	}
        pvs = lapply(subindices1, function(subinds){
          .processLine(splitline1,sums_1, subinds, splitline[1], collapse=collapse, 
                       debug=params$debug[[1]], lmTest=params$lmTest[[1]], canSkip=T)
        })
    
        if(FALSE){ ##this just to view the results for debugging
        
          df = lapply(pvs[ unlist( lapply(pvs, function(x)!is.null(x)))], function(x) x$pvs)
	       df1 = data.frame(t(data.frame(df)))
	      names(df1) = h_out[1:ncol(df1)]
        }
      },
      error=function(cond) {
        message(cond);
	break;
      }
    )
 
   # if(is.null(line)) break;
    #print(k)
  }
},
error=function(cond) {
  message(cond);
}
)
  closeAllConnections()
}



.checkMissing<-function(){
  miss1 = header[which(is.na(match(header, samples$NAME)))]
  t1 = table(unlist(lapply(miss1, function(x)strsplit(x,"-")[[1]][2])))
 print(t1)
  
  missing = samples[which(is.na(match(samples$NAME, header))),]
  missing$Cell_Type = factor(missing$Cell_Type)
  missing$Cohort = factor(missing$Cohort)
  missing$Cell_State = factor(missing$Cell_State)
  if(nrow(missing)>1){
    print("missing some barcodes..")
    print(table(missing$Cell_Type))
    print(table(missing$Cell_State))
    print(table(missing$Cohort))
  }
}

.readDefs<-function(defs_file){
  defs = read_json(defs_file) ## we split do we dont run into problems with too many connections
 # defs_all1 = lapply(defs_all, function(defs){
    for(i in 1:length(defs)) attr(defs[[i]], "nme")<-names(defs)[i]
    excl=unlist(lapply(defs, function(def) if(is.null(def$exclude)) FALSE else def$exclude))
  defs1 =   defs[!excl]
  #})
  
#  defs_all2 = lapply(defs_all1, function(defs){
      defs2 = unlist(lapply(defs1, function(def){
                .explode(def,levels(samples1[[def$type]]))
        }), rec=F)
      
 # })
  
#  defs_all2
      defs2
}

#explode(defs[[1]], levels(samples1[[defs[[1]]$type]]))
.explode<-function(def, levs){
  if(is.null(def$all_vs_all)) return(list(def))
  res = lapply(levs, function(lev){
    def1 = def
    def$case = lev
    def$control = NULL
    def
  })
  nme1 = attr(def,"nme")
  names(res )  = paste(nme1,levs,sep="_")
  for(i in 1:length(res)) attr(res[[i]], "nme")<-paste(nme1,levs[i],sep="_")
  res
}


.getSums<-function( input_file = getOption("input_file"), sum_file = getOption("sum_file","sums.rds")){
  #if(is.null(params$prefix))params$prefix=""
  closeAllConnections()
 # sum_file = paste0(getOption("prefix",""), "sums.rds")
  if(file.exists(sum_file)){
    print("reading sums")
    sums = readRDS(sum_file)
    return(sums)
    #sum_genes =   if(file.exists(sum_gene_file)) readRDS(sum_gene_file) else c()
  }else{
    connection = gzfile(input_file,"rb")
    # open(connection)
    line=readLines(connection,1)
    header=strsplit(line,split)[[1]]
    line = readLines(connection,1)
    sums = rep(0, length(header))
    sum_genes  = NULL
    i=0
    while(length(line)>0 && i< max_lines){
      row = strsplit(line,split)[[1]]
      
      if(i %% 1000 ==0){
        print(paste(row[1], i, sum(sums)))
      }
      if(length(row)==0) break;
      tryCatch(
        {
          
          
          inds=which(row!="0")[-1]
          if(length(inds)>0){
            rowv = as.numeric(row[inds])
            sums[inds]=sums[inds]+rowv
            
          }
        },
        error=function(cond) {
          print("ERROR")
          print(paste(length(inds), row[1], i))
          message(cond)
        }
      )
      line = readLines(connection,1)
      i= i+1
    }
    #print("SAVING IMAGE" )
    #save.image()
    max_lines=i+1
    saveRDS(sums, file = sum_file)
    # saveRDS(sum_genes, file = sum_gene_file)
  }
  sums[which(sums<getOption("back_thresh",0))] = NA
  names(sums) = header
  
  print(sum(sums))
  return(sums)
}


.sumTranscripts<-function( input_file, sum_file){
  if(file.exists(sum_file)) {
    print("sum file already exists")
    return();
  }
  split=getOption("split","\t")
  #if(is.null(params$prefix))params$prefix=""
  closeAllConnections()
  outf = gzfile(sum_file, open="w")
    connection = gzfile(input_file,"rb")
    # open(connection)
    line=readLines(connection,1)
    writeLines(line, outf)
    nxt = readLines2(connection,1,NULL) #
    k=0
    k1=0
    max_lines = getOption("max_lines",3e6) 
    while(!is.null(nxt) && k<max_lines){
      gene = nxt$gene
      nxt1 = readLines2(connection, 1, nxt)
      nxt=nxt1
      sums2 = .sumLines(nxt1$lines)
      all_sum = sum(sums2,na.rm=T)
      print(paste(k,gene,k1, all_sum))
      
      sums2[1] = gene
      line1=paste(sums2, collapse=split)
      writeLines(line1, outf)
      k = k+1
      k1=k1+length(nxt1$lines)
   }
  closeAllConnections() 
}


#NEW FUNCTS


.merge1<-function(t1,num_cols = c(),uniq_cols=c(), addName=NULL){
  
  t = t1[unlist(lapply(t1, function(x)!is.null(x)))]
  if(length(t)==0) return(NULL)
  t1 = t[[1]]
  if(!is.null(addName)){
    nme = rep(names(t)[[1]], nrow(t1))
    t1 = cbind(t1, nme)
  }
  if(length(t)>1){
    for(i in 2:length(t)){
      t_i=t[[i]]
      if(!is.null(addName)){
        nme = rep(names(t)[[i]], nrow(t_i))
        t_i = cbind(t_i, nme)
      }
      t1 = rbind(t1,t_i)
    }
  }
  dimnames(t1)[[1]]=1:nrow(t1)
  t2 = data.frame(t1)
  if(!is.null(addName)){
    names(t2)[ncol(t2)]=addName
  }
  
  if(length(uniq_cols)>0){
    facts = factor(apply(t2[,which(names(t2) %in% uniq_cols),drop=F],1,paste, collapse="."))
    t2 = t2[!duplicated(facts),,drop=F]
  }
  lev_cols =names(t2)[!( names(t2)%in% num_cols)]
  num_cols1 = num_cols[which(num_cols %in% names(t2))]
  for(j in num_cols1){
    
    t2[[j]] = as.numeric(t2[[j]])
  }
  for(j in lev_cols){
    t2[[j]] = factor(t2[[j]])
    
  }
  
  t2
}


.getCaseInds<-function(defcase,samples){
  case_inds =apply(data.frame( lapply( 1:length(defcase), function(j){ 
    if(defcase[[j]]=="all") return(rep(T, nrow(samples)))
    if(defcase[[j]]=="none") return(rep(F, nrow(samples)))
    samples[[names(defcase)[[j]]]] %in% defcase[[j]]
  } ))
  ,1,max)
}

.getDataFrame<-function(def,samples, header=NULL){
  case_inds = .getCaseInds(def$case, samples)
  control_inds =.getCaseInds(def$control,samples)
                    
  casecontrol = rep(NA, nrow(samples))
  casecontrol[case_inds==1]=1
  casecontrol[control_inds==1]=0
  inds1 = which(!is.na(casecontrol))
  casecontrol = casecontrol[inds1]
  index = 1:nrow(samples)
  if(is.null(header)){
    index_header = index0
  }else{
    index_header = match(samples$NAME, header)
  }
  index = index[inds1]
  index_header= index_header[inds1]
  all=rep("all",length(inds1))
  
  xT = rep(NA, length(inds1))
  xG = rep(NA, length(inds1))
  out= cbind(samples[inds1,,drop=F], index,index_header,casecontrol,all,xT,xG)
  out
}

.splitFrame<-function(splitBy ,df1 ){
  levs = unique(df1[[splitBy]]) 
  outd =lapply(levs, function(lev){
    df1[df1[[splitBy]]==lev,,drop=F]
  })
  names(outd)=levs
  
  outd
}

.splitFrameMultiple<-function(def, samples){
  dfs = list(samples)
  for(j in 1:length(def$splitBy)){
    dfs =  unlist(lapply(dfs, function(df1){
      splitBy=def$splitBy[[j]]
      .splitFrame(splitBy,df1)
    }),rec=F)
    for(k in 1:length(dfs)) attr(dfs[[k]],"nme")=names(dfs)[[k]]
  }
  dfs
}


.getFormulaFromDef<-function(def){
  str0=if(length(def$adjust_fixed)==0) "casecontrol~1" else paste0("casecontrol~",paste(unlist(c(def$adjust_fixed)),collapse="+"))
  str1=paste0("casecontrol~",paste(unlist(c("xG",def$adjust_fixed)),collapse="+"))
  str2=paste0("casecontrol~",paste(unlist(c("xT+xG",def$adjust_fixed)),collapse="+"))
  list(
    formula0=as.formula(str0),
    formula1=as.formula(str1),
    formula2 = as.formula(str2))
}
.pseudoBulk<-function(samples3,psb, def,
                      sumTags=c("Count"),
                      num_cols = c(sumTags,"index","index_header","xG","xT","casecontrol")){
  index_within=1:nrow(samples3)
  cell_type=attr(samples3,"nme")
  num_cols  = c(num_cols, "index_within")
  samples3 = cbind(samples3,index_within)
  s3 = factor( samples3[[psb]])
  levs = levels(s3)
  s4= lapply(levs, function(lev){
    inds = which(s3==lev)
    samples3[inds,,drop=F]
  })
  names(s4)=levs
  datf = data.frame(lapply(s4, function(s41){
    apply(s41, 2, function(x) length(unique(x)))
  }))
  dimnames(datf)[[1]] = names(s4[[1]])
  inds_to_keep=(apply(datf, 1, max)==1)
  sumInds =which( names(inds_to_keep) %in% sumTags)
  inds_to_keep[sumInds]=T
  inds_to_keep[which(names(inds_to_keep) %in% c("casecontrol"))]=TRUE
  s4_summ=
    .merge1(
      lapply(s4,function(s41){
        s41[1,sumInds] = apply(s41[,sumInds,drop=F],2,sum)
        s41[1,inds_to_keep]
      }), num_cols = num_cols)
  dimnames(s4_summ)[[1]] = names(s4)
  s4_sep = lapply(s4,function(s41){
    s41[,!inds_to_keep,drop=F]
  })
  inds2 = which(!(names(samples3) %in% num_cols))
  for(j in inds2) samples3[,j] = factor(samples3[,j])
  params_ = .getFormulaFromDef(def)
  gm0 = glm(params_$formula0, data=samples3,family="binomial")
  ll0 = logLik(gm0)

  gm0_1 = glm(params_$formula0, data=s4_summ,family="binomial")
  ll0_1 = logLik(gm0_1)  
  
  list(df_summ=s4_summ, sep=s4_sep, df_all=samples3,def=def,formula=params_,ll0=ll0,ll0_sum=ll0_1,
       cell_type=cell_type)
}


.addSplitColumns<-function(samples,cohortSplitName, cohortSplit="_"){
  if(!is.null(cohortSplit)){
    if(nchar(cohortSplit)>0){
      toSplit=names(cohortSplitName)[[1]]
      nmes1 = strsplit(cohortSplitName[[1]],":")[[1]]
      cohort1= lapply(as.character(samples[[toSplit]]), function(x){
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
  samples
}


.addMergeColumns<-function(samples, mergeColumns, cohortSplit="_"){
  if(!is.null(mergeColumns)){
    mergeColumns = getOption("mergeColumns")
    new_samples = data.frame(matrix(nrow = nrow(samples), ncol = length(mergeColumns)))
    names(new_samples) = names(mergeColumns)
    for(k in 1:length(mergeColumns)){
      new_samples[,k] = apply(samples[,match(mergeColumns[[k]], names(samples) )],1,paste,collapse=cohortSplit)
    }
    samples = cbind(samples, new_samples)
  }
  samples
}

.addClinicalData<-function(samples, clinical_data){
  if(!is.null(clinical_data)){
    clinical_data = read.table(clinical_data[[1]], sep="\t", header=T)
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
  samples
}


.addCountData<-function(samples){
  ##following makes a new meta file with the counts
  outf1 = paste0(params$meta_file,".1.txt")
  if(!file.exists(outf1)){
    sums = .getSums()
    names(sums) = header
    mi1 = match( samples$NAME,header)
    Count = sums[mi1]
    # indsNA1 = which(is.na(names(Count)))
    names(Count) = dimnames(samples)[[1]]
    samples1 = cbind(samples, Count)
    write.table(samples1, file=outf1,sep="\t",quote=F, row.names=F)
    samples = samples1
  }else{
    stop(paste(" redo using  ",outf1))
  }
  options("sumsInMeta"="Count")
  samples
}



.readMatrix<-function(connection, currLine1=NULL){
  geneColumn=getOption("geneColumn",2)
  # transcriptColumn =getOption("transcriptColumn",1)
  res = list()
  gene=NULL
  while(TRUE){
    if(is.null(currLine1)){
      break;
    }
    gene1=strsplit(currLine1[[1]],"\\|")[[1]][geneColumn]
    if(is.null(gene)) gene = gene1
    if(gene1!=gene) break;
    #  transcript = strsplit(currLine,"\\|")[[1]][transcriptColumn]
    which(currLine1!="0")
    res[[length(res)+1]] = currLine1
    currLine=readLines(connection,1)
    currLine1=if(is.null(currLine)) NULL else strsplit(currLine,params$split[[1]])[[1]]
  }
  names(res) = lapply(res, function(r)r[[1]])
  attr(res,"next")=currLine1
  attr(res,"gene")=gene
  res
}



.assocTestAll<-function(psb1, matr,
                        min_num_reads = getOption("min_num_reads",20)
){
  df_all = psb1$df_all
  df_summ=psb1$df_summ
  gene=attr(matr,"gene")
  matr1 = data.frame(lapply(matr, function (line) as.numeric(line[psb1$df_all$index_header])))
  matr2 =  t(data.frame(lapply(psb1$sep, function(sep1){
    apply(matr1[sep1$index_within,,drop=F],2,sum)
  })))
  
  df_all$xG = apply(matr1,1,sum)
  df_summ$xG = apply(matr2,1,sum)
  tsums=apply(matr2,2,sum) ## can use matr2 or matr1
  inds_t = which(tsums>min_num_reads)
  matr1 = matr1[,inds_t,drop=F]
  matr2 = matr2[,inds_t,drop=F]
  pv_res=.sigInner(df_all, matr1,params_, psb1$ll0,"cell")
  pv_res2=.sigInner(df_summ, matr2,params_,psb1$ll0_sum,"pb")
  Name = factor(rep(attr(psb1$def,"nme"), nrow(pv_res)))
  Condition = paste(paste(def$case, collapse="_"), "vs",paste(def$control, collapse="_"),sep="_")
  Cell_type = factor(rep(psb1$cell_type, nrow(pv_res)))
  ID = factor(dimnames(pv_res)[[1]])
  result=cbind(ID, Name,Condition,Cell_type, pv_res, pv_res2)
  result
}

.sigInner<-function(data,matr, params_, ll0, type=""){
  ##SIGNIFICANCE AT CELL LEVEL
  formula1=params_$formula1; formula2=params_$formula2
  gm1 = glm(formula1, data=data,family="binomial")
  ll1 = logLik(gm1)
  pv_res = matrix(.chisq(ll0,gm1),nrow=1)  
  dimnames(pv_res) =list(gene,paste(c("pv","beta"),type,sep="_" ))#,"count","non-zero") 
  if(nrow(matr)==0) return(pv_res)
  pv_res1 = t(apply(matr,2,function(v){
    data$xT = v
    gm2 = glm(formula2, data=data,family="binomial")
    .chisq(ll1,gm2)
  }))
  pv_res2 = data.frame(rbind(pv_res,pv_res1))
  pv_res2
}


.chisq<-function(ll1,gm2){
  s1 = summary(gm2)
  ll2 =logLik(gm2)
  df = abs(attr(ll1,"df")[1]-attr(ll2,"df")[1])
  lldif= abs(2*(ll1 - ll2))
  pv = pchisq(lldif,df,lower.tail=FALSE,log.p=F)
  res = c(pv, s1$coefficients[2,1])
  res
}

