#' @title
#' oxBScut
#'
#' @name
#' oxBScut
#'
#' @description
#' The function makes a recommendation for a 5hmC cut-off to eliminate unstable
#' 5hmC probes on the Illumina HumanMethylation BeadChip and evaluate the 
#' information loss with the low-quality probes.
#'
#' @param
#' matrix_5hmc is a matrix or dataframe of 5hmC beta values.
#' poobha_5hmc is a matrix or dataframe of p-values with out-of-band array 
#' hybridization that is corresponding to the 5hmC beta matrix.
#'
#' @return
#' Recommended 5hmC cut-off, low-quality fraction, a figure on 5hmC distribution
#' between high-quality probe and low-quality probe, a figure on low-quality
#' CpG genomic location enrichment, a figure on coefficient of variation and 
#' median 5hmC with the recommended 5hmC cut-off marked.  
#'
#' @import dplyr
#' @import matrixStats
#' @import RnBeads
#' @import EnvStats
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#' @import epitools
#' @import forestplot
#' @import ggplot2
#' 
#' @export


oxBScut<-function (matrix_5hmc,poobha_5hmc) {
  library(dplyr)
  library(matrixStats)
  library(RnBeads)
  library(EnvStats)
  library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  median_5hmc<-as.data.frame(rowMedians(as.matrix(matrix_5hmc)))
  rownames(median_5hmc)<-rownames(matrix_5hmc)
  
  cv_5hmc<-as.data.frame(apply(matrix_5hmc, 1, cv)) 
  colnames(cv_5hmc)[1]<-"5hmC"
  
  cv_median_5hmc<-cbind(median_5hmc,cv_5hmc)
  colnames(cv_median_5hmc)<-c("MedianBeta","CV")
  
  
  cv_median_5hmc<-cv_median_5hmc[order(cv_median_5hmc$MedianBeta, decreasing = T),]
  ID<-intersect(rownames(poobha_5hmc),rownames(cv_median_5hmc))
  poobha_5hmc_1<-poobha_5hmc[rowSums(poobha_5hmc<0.05)>dim(poobha_5hmc)[2]*0.2,]
  cv_median_5hmc$LT<-(cv_median_5hmc$CV)^2 + (cv_median_5hmc$MedianBeta)^2
  cutoff<-cv_median_5hmc[cv_median_5hmc$LT == min(cv_median_5hmc$LT,na.rm = TRUE),"MedianBeta"][[1]]
  ID1<-rownames(cv_median_5hmc[cv_median_5hmc$MedianBeta>cutoff,])
  cv_median_5hmc$MBCut<-ifelse(rownames(cv_median_5hmc)%in%ID1, "Yes", "No")
  cv_median_5hmc$pCut<-ifelse(rownames(cv_median_5hmc)%in%rownames(poobha_5hmc_1),"Yes","No")
  c_Table <- as.data.frame(matrix(ncol = 3,nrow = 3))
  colnames(c_Table) <- c('#pOOBHA<0.05 >20%','#pOOBHA<0.05 <20%','Total')
  rownames(c_Table) <- c('Meidan Beta > cutoff','Median Beta < cutoff','Total')
  c_Table[1,1]<-table(cv_median_5hmc$MBCut,cv_median_5hmc$pCut)[2,2]
  c_Table[1,2]<-table(cv_median_5hmc$MBCut,cv_median_5hmc$pCut)[2,1]
  c_Table[2,1]<-table(cv_median_5hmc$MBCut,cv_median_5hmc$pCut)[1,2]
  c_Table[2,2]<-table(cv_median_5hmc$MBCut,cv_median_5hmc$pCut)[1,1]
  c_Table[3,1]<-sum(c_Table[1:2,1])
  c_Table[3,2]<-sum(c_Table[1:2,2])
  c_Table[1,3]<-sum(c_Table[1,1:2])
  c_Table[2,3]<-sum(c_Table[2,1:2])
  c_Table[3,3]<-sum(c_Table[3,1:2])
  LQP<-round(c_Table[1,2]/c_Table[1,3]*100,2)
  Removed_CpGs<-rownames(cv_median_5hmc[cv_median_5hmc$MedianBeta>cutoff&cv_median_5hmc$pCut=="No",])
  Background_cpGs<-ID1
  
  if (c_Table[3,3]>800000){
    annot<-getAnnotation("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
  } else {
  annot<-getAnnotation("IlluminaHumanMethylation450kanno.ilmn12.hg19")}
  
  Annotation <- annot[annot$Name %in% Background_cpGs,]
  
  Annotation$Relation_to_Island <- as.character(Annotation$Relation_to_Island)
  Annotation$islandAdjust <- ifelse(Annotation$Relation_to_Island %in% c('N_Shelf','S_Shelf'),'Shelf',
                                    ifelse(Annotation$Relation_to_Island %in% c('N_Shore','S_Shore'),'Shore',ifelse(Annotation$Relation_to_Island == "Island","Island","OpenSea")))
  
  Annotation$DMR <- ifelse(Annotation$Name %in% Removed_CpGs,1,0)
  
  #### Relation to CpG Island Enrichment Test
  #### Relation to CpG Island Enrichment Test
  # Calculate odds ratios
  orTable1 <- as.data.frame(matrix(ncol = 5,nrow = 0))
  colnames(orTable1) <- c('proportion','OR','lower','upper','pVal')
  islandLevels <- c('Island','Shore','Shelf','OpenSea')
  islandLabels <- c('Island','Shore','Shelf','OpenSea')
  library(epitools)
  for (i in islandLevels) {
    # Calculate proportion of sites in that context
    prop <- table(Annotation$DMR, Annotation$islandAdjust  == i)['1','TRUE']/sum(Annotation$DMR)
    # Calculate odds ratio
    tempOR <- oddsratio.fisher(table(Annotation$DMR,Annotation$islandAdjust  == i))
    # Add values to output
    orTable1[i,] <- c(prop,tempOR$measure[2,],tempOR$p.value[2,'fisher.exact'])
  }
  orText <- cbind(c('',islandLevels),
                  c('Odds Ratio (95% CI)',paste(format(round(orTable1$OR,3),nsmall = 2),' (',format(round(orTable1$lower,3),nsmall = 2),', ',
                                                format(round(orTable1$upper,3),nsmall = 2),')',sep = '')),
                  c('p value',formatC(orTable1$pVal,digits = 3)),
                  c('Proportion of Sites',format(round(orTable1$proportion,3),nsmall = 2)))
  
  library(forestplot)
  f1<-forestplot(labeltext = orText, graph.pos = 2, 
             mean = c(NA,orTable1$OR),
             lower = c(NA,orTable1$lower),
             upper = c(NA,orTable1$upper),
             title="Low Quality CpG Genomic Location Enrichment",
             
             txt_gp=fpTxtGp(label=gpar(cex=1),
                            ticks=gpar(cex=1),
                            xlab=gpar(cex = 1.9),
                            title=gpar(cex = 1.5)),
             
             #xticks = c(0,1,2),
             #xlog = T,
             col=fpColors(box="black", lines="black", zero = "gray50")
             ,
             fn.ci_norm = fpDrawPointCI, pch = 16
             , cex=1,
             zero = 1
             , 
             lineheight = "auto"
             ,
             boxsize=0.08
             , 
             colgap=unit(6,"mm")
  )
  
  removed<-stack(matrix_5hmc[Removed_CpGs,])
  removed$probe_type<-"Low Quality Probes"
  keep<-rownames(cv_median_5hmc[cv_median_5hmc$MedianBeta>cutoff&cv_median_5hmc$pCut=="Yes",])
  keep_5hmc<-stack(matrix_5hmc[keep,])
  keep_5hmc$probe_type<-"High Quality Probes"
  stack<-rbind(removed,keep_5hmc)
  library(ggplot2)
  f2<-ggplot(stack, aes(x=values,fill = probe_type, color=probe_type)) +
    geom_histogram(position="identity", alpha=0.5, bins = 50)+xlab("5hmC Beta")+ylab("Frequency")+
    #scale_x_continuous(breaks=seq(0, 0.8, 0.1), limits=c(0, 0.8))+
    theme(legend.position="top", axis.title = element_text( size = 15 , face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  f3<-create.densityScatter(cv_median_5hmc[,c("MedianBeta","CV")], sparse.points = 0)+xlab("5hmC Median Beta")+ylab("Coefficient of Variation")+
    #scale_y_continuous(limits=c(0, 4.3),breaks = seq(0, 4, by = 1))+scale_x_continuous(breaks = seq(0, 0.6, by = 0.1))+
    geom_vline(xintercept = cutoff, linetype="dotted", color = "red")+
    theme(axis.text = element_text( size = 10 ),
          axis.title = element_text( size = 15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  print(paste0("Recommended 5hmC cut-off is ",round(cutoff,5)))
  print(paste0("Low-quality fraction is ",LQP))
  print(f1)
  print(f2)
  print(f3)
}

