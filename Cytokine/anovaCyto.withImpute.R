require(plyr)
library(dplyr)
library(tidyverse)
library(broom)
library(ggplot2)

cyto <- read.table("cytokine_average_ANOVA_input.tsv", sep="\t",header=TRUE,row.name=1,na.strings="NA")

groups <- read.csv("./groups.csv",header=TRUE)
rownames(groups) <- groups$Sample
groups$Treatment <- factor(groups$Treatment, levels = c("PBS_XXM", "PBS_XXF", "PBS_XYM", "PBS_XYF", "HDM_XXM", "HDM_XXF", "HDM_XYM","HDM_XYF"))

tube2sample <- read.table("tube2sample.tsv", sep="\t",header=TRUE,row.name=1)
cyto2 <- merge(tube2sample,cyto,by=0,sort=F)
rownames(cyto2) <- cyto2$Sample
cyto2 <- subset(cyto2,select=-c(Row.names,Sample))


cyto_names <- names(cyto2)
cytoplot <- merge(cyto2,groups,by=0,sort=F)
cytoplot$TestGroup <- gsub("_Group_\\d+","",cytoplot$Experimental)
cytoplot$TestGroup <- factor(cytoplot$TestGroup, levels=c("Control","Test"))

logCyto <- log2(cyto2)
logcytoplot <- merge(logCyto,groups,by=0,sort=F)
logcytoplot$TestGroup <- gsub("_Group_\\d+","",logcytoplot$Experimental)
logcytoplot$TestGroup <- factor(logcytoplot$TestGroup, levels=c("Control","Test"))

anovaTests <- list()
cairo_pdf("cytokine_anova_plots.withImpute.pdf", onefile=TRUE)

for(i in cyto_names){
  print(i)
  ploti <- gsub("\\.","-",i)
  natest <- logcytoplot  %>% group_by(Treatment) %>% dplyr::summarize(count_na = sum(! is.na(!!sym(i))))
# Test is true if at least one sample is not NA
  if(sum(is.na(logcytoplot[[i]])) < 24 ){
#   Impute value minimum cytokine concentration
#   Note all Out of Range values were due to low fluorescence, so we can take the minimum 
    logcytoplot[[i]][is.na(logcytoplot[[i]])] <- log2(min(cytoplot[[i]],na.rm=TRUE))
    cytoplot[[i]][is.na(cytoplot[[i]])] <- min(cytoplot[[i]],na.rm=TRUE)
    print(paste("testing ",i,sep=""))
    test <- aov(logcytoplot[,i] ~ Treatment, data = logcytoplot)
    if(exists("dfFinal")){
      dfFinal <- rbind.fill(dfFinal,data.frame(cyto=i,broom::tidy(test)))
    }else{
      dfFinal <- data.frame(cyto=i,broom::tidy(test))
    }
#   test normality
    stest <- shapiro.test(test$residuals)
    if(exists("shapiro")){
      shapiro <- rbind.fill(shapiro,data.frame(cyto=i,broom::tidy(stest)))
    }else{
      shapiro <- data.frame(cyto=i,broom::tidy(stest))
    }
    par(mfrow=c(2,2))
    plot(test,main=ploti)
    par(mfrow=c(1,1))
    hist(test$residuals,main=paste(ploti," ANOVA Residuals Histogram",sep=""), xlab=paste(ploti," ANOVA Residuals",sep=""))
    anovaTests[[i]] <- test
  }
  print(ggplot(cytoplot, aes(x=Treatment,y=!!sym(i))) + geom_boxplot(aes(color=Treatment, fill=Treatment), alpha=0.3) + geom_point() + facet_wrap(vars(TestGroup),scales="free_x") + guides(alpha="none") + ylab("Observed Concentration (pg/mL)") + labs(title=ploti) + theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(size=8)))
  if(sum(is.na(logcytoplot[[i]])) > 23 ){
    logcytoplot[[i]] <- as.logical(logcytoplot[[i]])
  }
  print(ggplot(logcytoplot, aes(x=Treatment,y=!!sym(i))) + geom_boxplot(aes(color=Treatment, fill=Treatment), alpha=0.3) + geom_point() + facet_wrap(vars(TestGroup),scales="free_x") + guides(alpha="none") + ylab("log2(Observed Concentration (pg/mL))") + labs(title=paste("log2(",ploti,")",sep="")) + theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(size=8)))
}
write.table(dfFinal,file="cytokineANOVA.impute.tsv", quote=FALSE, sep="\t", na="NA", row.names=FALSE)
write.table(shapiro,file="cytokineANOVA.shapiro.impute.tsv", quote=FALSE, sep="\t", na="NA", row.names=FALSE)
dev.off()
save(anovaTests,cytoplot,logcytoplot,file="anovaTests.impute.RData")
writeLines(capture.output(sessionInfo()), "sessionInfo.cytokineANOVA.impute.txt")
