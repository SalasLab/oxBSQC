#' @title
#' oxBSlabel
#'
#' @name
#' oxBSlabel
#'
#' @description
#' The function predicts the BS or oxBS label for paried bisulfite and oxidative
#' bisulfite treated samples using the internal normalization control probes on
#' the Illumina HumanMethylation BeadChip. 
#' 
#'
#' @param
#' rgset 'RGChannelSet' including pheno matrix with columns "Sample_ID" that is
#' unique to each sample and “Subject" that is unique to each paired subject.
#'
#'
#'
#' @return
#' A pheno matrix with columns on predicted BS and oxBS label.  
#'
#' @import 	minfi
#'
#'
#' @export


oxBSlabel<-function (rgSet) {
  if ( !is(rgSet, "RGChannelSet")) {
    stop("object needs to be of class 'RGChannelSet'")
  }


    ctrls <- getProbeInfo(rgSet, type = "Control")
  
  ctrls <- ctrls[ctrls$Address %in% rownames(rgSet), ]
  ctrl_r <- assays(rgSet)$Red[ctrls$Address, ]
  ctrl_g <- assays(rgSet)$Green[ctrls$Address, ]
  cc <- ctrls[ctrls$Type %in% c("NORM_A", "NORM_C",
                                "NORM_G", "NORM_T"), ]
  red <- ctrl_r[cc$Address, ]
  grn <- ctrl_g[cc$Address, ]
  ymax <- max(red, grn) * 1.01
  
  ccdf<-as.data.frame(cc)
  #identical(rownames(grn),ccdf$Address)
  NORM_A<-ccdf %>% filter(Type == "NORM_A")
  NORM_G<-ccdf %>% filter(Type == "NORM_G")
  NORM_C<-ccdf %>% filter(Type == "NORM_C")
  NORM_T<-ccdf %>% filter(Type == "NORM_T")
  
  NORM_A_Red<-red[NORM_A$Address,]
  NORM_C_Green<-grn[NORM_C$Address,]
  NORM_G_Green<-grn[NORM_G$Address,]
  NORM_T_Red<-red[NORM_T$Address,]
  
  pheno<-as.data.frame(rgSet@colData)
  pheno<-pheno[order(pheno$Subject),]
  
  NORM_A_Red<-NORM_A_Red[,pheno$Sample_ID]
  NORM_C_Green<-NORM_C_Green[,pheno$Sample_ID]
  NORM_G_Green<-NORM_G_Green[,pheno$Sample_ID]
  NORM_T_Red<-NORM_T_Red[,pheno$Sample_ID]
  
  P_Matrix<-matrix(NA,nrow = nrow(pheno)/2,ncol = 5)
  
  rownames(P_Matrix)<-unique(pheno$Subject)
  colnames(P_Matrix)<-c("ID","Norm_A Red P-value","Norm_C Green P-value",
                        "Norm_G Green P-value","Norm_T Red P-value")
  
  even <- seq(2,nrow(pheno),2)
  dd <- seq(1,by=2, len=nrow(pheno)/2)
  
  for (i in 1:(nrow(pheno)/2)) {
    
    P_Matrix[i,2]<-wilcox.test(NORM_A_Red[,dd[i]],NORM_A_Red[,even[i]], alternative = "greater", paired = TRUE)$p.value
    
  }
  
  for (i in 1:(nrow(pheno)/2)) {
    
    P_Matrix[i,3]<-wilcox.test(NORM_C_Green[,dd[i]],NORM_C_Green[,even[i]], alternative = "greater", paired = TRUE)$p.value
    
  }
  
  for (i in 1:(nrow(pheno)/2)) {
    
    P_Matrix[i,4]<-wilcox.test(NORM_G_Green[,dd[i]],NORM_G_Green[,even[i]], alternative = "greater", paired = TRUE)$p.value
    
  }
  
  for (i in 1:(nrow(pheno)/2)) {
    
    P_Matrix[i,5]<-wilcox.test(NORM_T_Red[,dd[i]],NORM_T_Red[,even[i]], alternative = "greater", paired = TRUE)$p.value
    
  }
  
  P_Matrix[,1]<-rownames(P_Matrix)
  P_Matrix<-as.data.frame(P_Matrix)
  
  for (i in 2:5) {
    P_Matrix[,i]<-as.numeric(P_Matrix[,i])
  }
  
  
  c<-c()
  
  for (j in 1:(nrow(pheno)/2)) {
    
    
    x<-ifelse(P_Matrix[j,2]<0.05, "BS","oxBS" )
    c<-c(c,x)
  }
  
  pheno[dd,ncol(pheno)+1]<-c  
  
  for (i in even) {
    pheno[i,ncol(pheno)]<-ifelse(pheno[i-1,ncol(pheno)]=="oxBS","BS","oxBS")
  }
  #----------------------------------------------------------------------------
  c<-c()
  
  for (j in 1:(nrow(pheno)/2)) {
    
    
    x<-ifelse(P_Matrix[j,3]<0.05, "BS","oxBS" )
    c<-c(c,x)
  }
  
  pheno[dd,ncol(pheno)+1]<-c  
  
  for (i in even) {
    pheno[i,ncol(pheno)]<-ifelse(pheno[i-1,ncol(pheno)]=="oxBS","BS","oxBS")
  }
  #-----------------------------------------------------------------------------
  c<-c()
  
  for (j in 1:(nrow(pheno)/2)) {
    
    
    x<-ifelse(P_Matrix[j,4]<0.05, "BS","oxBS" )
    c<-c(c,x)
  }
  
  pheno[dd,ncol(pheno)+1]<-c  
  
  for (i in even) {
    pheno[i,ncol(pheno)]<-ifelse(pheno[i-1,ncol(pheno)]=="oxBS","BS","oxBS")
  }
  #----------------------------------------------------------------------------
  c<-c()
  
  for (j in 1:(nrow(pheno)/2)) {
    
    
    x<-ifelse(P_Matrix[j,5]<0.05, "BS","oxBS" )
    c<-c(c,x)
  }
  
  pheno[dd,ncol(pheno)+1]<-c  
  
  for (i in even) {
    pheno[i,ncol(pheno)]<-ifelse(pheno[i-1,ncol(pheno)]=="oxBS","BS","oxBS")
  }
  
  colnames(pheno)[(ncol(pheno)-3):ncol(pheno)]<-c("Norm_A Red Predicted","Norm_C Green Predicted","Norm_G Green Predicted","Norm_T Red Predicted")
  
  return(pheno)
  
}


