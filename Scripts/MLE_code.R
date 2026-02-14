library(dplyr)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(parallel)

### tidy the cellsnp result
fun.loci2pare <-function(Parents, Sample){ 
  if(Parents=="ULd2"){
    parent.chrom.vcfpath <-c("forR/UL.parents.all.forR")
  } else if(Parents=="LUd2"){
    parent.chrom.vcfpath <-c("forR/LUd2.parents.all.forR")
  }else{
    parent.chrom.vcfpath <-c("forR/LUd7.parents.all.forR")
  }
  parent.chrom.VCF <-read.table(parent.chrom.vcfpath, 
                                col.names=c("Posi","TL_base","TU_base"))
  
  chrom.cellsnp.path <-paste0(Sample,".cell2snp.forR")
  chrom.cellsnp <- read.table(chrom.cellsnp.path, 
                              col.names =c("CHROM","POS","Posi","Ref","Alt","Ref_count","Alt_count","GT","Barcode"))
  
  chrom.Count <- inner_join(chrom.cellsnp, parent.chrom.VCF, by="Posi")

  tmp11 <-
    chrom.Count %>% 
    mutate(TL_count = if_else(TL_base == Ref, Ref_count, Alt_count), TU_count = if_else(TU_base == Ref, Ref_count, Alt_count)) 
  
  var_name <- paste(Sample, "cell2SNP", sep = ".")  
  assign(var_name, tmp11, envir = .GlobalEnv)  
}

fun.loci2pare("ULd2","UL1")
fun.loci2pare("ULd2","UL2")
fun.loci2pare("LUd2","LU1")
fun.loci2pare("LUd2","LU2")
fun.loci2pare("LUd7","LU3")
fun.loci2pare("LUd7","LU4")

UL1.SNP.forInte <- UL1.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_5")) %>% mutate(orig.ident="UL1") 
UL2.SNP.forInte <- UL2.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_6")) %>% mutate(orig.ident="UL2") 
LU1.SNP.forInte <- LU1.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_3")) %>% mutate(orig.ident="LU1")
LU2.SNP.forInte <- LU2.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_4")) %>% mutate(orig.ident="LU2") 
LU3.SNP.forInte <- LU3.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_1")) %>% mutate(orig.ident="LU3") 
LU4.SNP.forInte <- LU4.cell2SNP %>% mutate(Barcode=paste0(Barcode,"_2")) %>% mutate(orig.ident="LU4") 
inte.SNP <- rbind(UL1.SNP.forInte, UL2.SNP.forInte, LU1.SNP.forInte, LU2.SNP.forInte, LU3.SNP.forInte, LU4.SNP.forInte)

inte.U5 <- inte.SNP %>% inner_join(inte.obj@meta.data %>% select(Barcode,tissue),by='Barcode') %>% select(Posi, orig.ident, Barcode, tissue, TL_base, TU_base, TL_count, TU_count) %>% 
  filter(TL_count + TU_count >= 5)

UL1.U5C5 <- inte.U5 %>% filter(orig.ident=="UL1") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()
UL2.U5C5 <- inte.U5 %>% filter(orig.ident=="UL2") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()
LU1.U5C5 <- inte.U5 %>% filter(orig.ident=="LU1") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()
LU2.U5C5 <- inte.U5 %>% filter(orig.ident=="LU2") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()
LU3.U5C5 <- inte.U5 %>% filter(orig.ident=="LU3") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()
LU4.U5C5 <- inte.U5 %>% filter(orig.ident=="LU4") %>% group_by(Posi) %>% filter(n() >= 5) %>% ungroup()

UL1.splitdata <- UL1.U5C5 %>% split(x=.,f=.$Posi)
UL2.splitdata <- UL2.U5C5 %>% split(x=.,f=.$Posi)
LU1.splitdata <- LU1.U5C5 %>% split(x=.,f=.$Posi)
LU2.splitdata <- LU2.U5C5 %>% split(x=.,f=.$Posi)
LU3.splitdata <- LU3.U5C5 %>% split(x=.,f=.$Posi)
LU4.splitdata <- LU4.U5C5 %>% split(x=.,f=.$Posi)

### MLE
LH.Model_0 <- function(DDData){
  tmp.LH <- 
    DDData %>% 
    mutate(logP=dbinom(Read1,Read1+Read2,0.5) %>% log) %>% 
    .$logP %>% sum
  
  return(tmp.LH)
}

LH.Model_1 <- function(DDData,p){
  tmp.LH <- 
    DDData %>% 
    mutate(logP=dbinom(Read1,Read1+Read2,p) %>% log) %>% 
    .$logP %>% sum
  
  return(tmp.LH)
}

LH.Model_2 <- function(DDData,p,d_vector){
  tmp.LH <- 
    DDData %>% 
    mutate(logP=dbinom(Read1,Read1+Read2,if_else(d_vector==1,p,1-p)) %>% log) %>% 
    .$logP %>% sum
  
  return(tmp.LH)
}

Fun.LRT <- function(LH_1,LH_2,DF){
  # Compute the LRT statistic
  LRT_stat <- 2 * (LH_2 - LH_1)
  
  # Compute the p-value from the chi-squared distribution
  p_value <- 1-pchisq(LRT_stat, df = DF)
  
  return(p_value)
}


Fun.RunMLE <- function(DDData){
  tmp1 <- DDData
  
  ## Estimation
  # Model 0, f=0
  message("Estimating Model 0")
  res.model0 <- tmp1 %>% LH.Model_0
  
  # Model 1, f=1
  message("Estimating Model 1")
  tmp.initiate_1 <- tmp1 %>% mutate(PPP=Read1/(Read1+Read2)) %>% .$PPP %>% mean
  res.model1 <- 
    optim(tmp.initiate_1,fn=function(ppp){ 
      tmp.out <- -LH.Model_1(tmp1,ppp) 
      if(is.infinite(tmp.out)){ tmp.out <- 99999 }
      return(tmp.out)
    },method="L-BFGS-B",lower = c(1e-5),upper = c(1-1e-5))
  
  # Model 2, f=2
  message("Estimating Model 2")
  
  
  tmp.set <- 
    list(
      p=seq(1e-5,1-1e-5,0.1),
      q=seq(1e-5,1-1e-5,0.1)
    ) %>% 
    expand.grid
  
  tmp.initiate_2 <- 
    1:nrow(tmp.set) %>% 
    lapply(function(iii){
      tmp.d_vector <- tmp1 %>% mutate(d=if_else(percent_rank(Read1/(Read1+Read2))>tmp.set$q[iii],1,-1)) %>% .$d
      tmp.LH <- tmp1 %>% LH.Model_2(p=tmp.set$p[iii],d_vector = tmp.d_vector)
      
      data.frame(
        p=tmp.set$p[iii],
        q=tmp.set$q[iii],
        LH=tmp.LH
      )
    }) %>% bind_rows %>% 
    filter(min_rank(-jitter(LH))==1) %>% 
    select(p,q) %>% as.numeric    
  
  res.model2 <- 
    optim(tmp.initiate_2,fn=function(VVVector){ 
      ppp <- VVVector[1]
      qqq <- VVVector[2]
      tmp.d_vector <- tmp1 %>% mutate(d=if_else(percent_rank(Read1/(Read1+Read2))>qqq,1,-1)) %>% .$d
      tmp.out <- -LH.Model_2(tmp1,p=ppp,d_vector = tmp.d_vector)
      if(is.infinite(tmp.out)){ tmp.out <- 99999 }
      return(tmp.out)
    },method="L-BFGS-B",lower = c(1e-5,1e-5),upper = c(1-1e-5,1-1e-5))
  
  ## Maximum likelihood
  message("Calculating Maximum Likelihood")
  tmp.estimation <- 
    list(
      M0=data.frame(Model="M0",p=0.5,q=1,LH=res.model0),
      M1=data.frame(Model="M1",p=res.model1$par,q=1,LH=-res.model1$value),
      M2=data.frame(Model="M2",p=res.model1$par[1],q=res.model1$par[2],LH=-res.model2$value)
    )
  
  ## Likelihood ratio test
  message("Performing Likelihood Ratio Tests")
  tmp.LRT <- 
    list(
      M0_M1=
        Fun.LRT(
          res.model0,
          -res.model1$value,
          1
        ),
      M1_M2=
        Fun.LRT(
          -res.model1$value,
          -res.model2$value,
          1
        )
    )
  
  tmp.bestmodel <- if_else(tmp.LRT$M0_M1<0.05,if_else(tmp.LRT$M1_M2<0.05,"M2","M1"),"M0")
  
  ## Return
  tmp.out <- list(
    BestModel=tmp.bestmodel,
    Estimates=tmp.estimation[[tmp.bestmodel]],
    LRT=tmp.LRT,
    History=list(
      M0=res.model0,
      M1=res.model1,
      M2=res.model2
    )
  )
  
  return(tmp.out)
}


### LU1
tmp.dropLU1 <- 
  names(LU1.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/LU1/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/LU1/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        LU1.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })

### LU2
tmp.dropLU2 <- 
  names(LU2.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/LU2/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/LU2/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        LU2.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })


### LU3
tmp.dropLU3 <- 
  names(LU3.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/LU3/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/LU3/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        LU3.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })

### LU4
tmp.dropLU4 <- 
  names(LU4.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/LU4/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/LU4/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        LU4.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })

### UL1
tmp.dropUL1 <- 
  names(UL1.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/UL1/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/UL1/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        UL1.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })

### UL2
tmp.dropUL2 <- 
  names(UL2.splitdata) %>% #head(1) %>% 
  mclapply(mc.cores = 20,function(mmm){
    # lapply(function(mmm){
    message(mmm)
    
    tmp.chr <- str_split_1(mmm,"_") %>% .[[1]]
    system(paste0("mkdir -p MLE_result_byIndiv/UL2/",tmp.chr))
    
    tmp.fileout <- paste0("MLE_result_byIndiv/UL2/",tmp.chr,"/",mmm,".RData")
    
    if(!file.exists(tmp.fileout)){
      tmp1 <- 
        UL2.splitdata[[mmm]] %>% 
        group_by %>% 
        select(Barcode,Read1=TL_count,Read2=TU_count) # using m/p or TL/TU
      
      tmp.res <- Fun.RunMLE(tmp1)
      
      save(tmp.res,file=tmp.fileout)
    }else{
      message("Skip")
    }
    
    return()
  })


### P adjust
Func_padjust <- function(Sample_Name, data_dir = "MLE_result_byIndiv/", cores = 10) {
  tmp.filelist <- dir(file.path(data_dir, Sample_Name), recursive = TRUE, full.names = TRUE)
  
  tmp.multitest <-
    tmp.filelist %>%
    mclapply(mc.cores = cores,function(fff){
      load(fff)
      tmp.posi <-  str_extract(fff, "[^/]+\\.RData") %>% str_remove("\\.RData")
      data.frame(Posi = tmp.posi,
                 M0_M1_p = tmp.res$LRT$M0_M1,
                 M1_M2_p = tmp.res$LRT$M1_M2) 
    }) %>% bind_rows()
  
  tmp1 <- 
    tmp.multitest %>% 
    mutate(M0_M1_padj=p.adjust(M0_M1_p,"bonferroni"),
           M0_M1_sig=M0_M1_padj<0.05)
  
  tmp2 <- 
    tmp1 %>% filter(M0_M1_sig) %>% 
    mutate(M1_M2_padj=p.adjust(M1_M2_p,"bonferroni"),
           M1_M2_sig=M1_M2_padj<0.05)
  
  tmp.adjmodel <- 
    rbind(
      tmp1 %>% filter(!M0_M1_sig) %>% mutate(Model="M0") %>% select(Posi,Model),
      tmp2 %>% mutate(Model=if_else(M1_M2_sig,"M2","M1")) %>% select(Posi,Model)
    )
  
  
  tmp.MLEdata <- 
    tmp.filelist %>%
    mclapply(mc.cores = cores,function(fff){
      load(fff)
      tmp.posi <-  str_extract(fff, "[^/]+\\.RData") %>% str_remove("\\.RData")
      
      tmp.res$Estimates %>% 
        mutate(Posi=tmp.posi) %>% 
        group_by(Posi,Model,p,LH) %>% 
        summarise(q=mean(q),n=n(),
                  Ncells=get(paste0(Sample_Name, ".splitdata"))[[tmp.posi]] %>% nrow(),.groups = "drop") %>%
        select(Posi,Ncells,Model,p,q,n,LH)
    }) %>% bind_rows()
  
  Modelinfo <- tmp.adjmodel %>% 
    left_join(tmp.MLEdata %>% dplyr::select(-Model), by = 'Posi') %>% 
    mutate(Indiv_rename = Sample_Name)
  
  var_name <- paste(Sample_Name, "Modelinfo", sep = ".")  
  assign(var_name, Modelinfo, envir = .GlobalEnv)  
  
}

Func_padjust("UL1")
Func_padjust("UL2")
Func_padjust("LU1")
Func_padjust("LU2")
Func_padjust("LU3")
Func_padjust("LU4")
