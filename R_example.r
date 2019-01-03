library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)


pheno = read.csv("/home/joshij/Desktop/Sample_II/DATA_table_phenotype.csv")

bg = ballgown(dataDir='/home/joshij/Desktop/Sample_II/Result_files/Merged_STRINGTIE_Out_files', samplePattern = 'SRR', pData=pheno)

bg_filt  = subset(bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)


results_trans = stattest(bg_filt, feature="transcript", covariate="SEX", adjustvars=c("Organs"), getFC=TRUE, meas="FPKM")

results_genes = stattest(bg_filt, feature="gene", covariate="SEX", adjustvars=c("Organs"), getFC=TRUE, meas="FPKM")

results_transcript = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),results_trans)

results_transcript = arrange(results_transcript,pval)

results_genes = arrange(results_genes,pval)

write.csv(data_frame, file_name, now.names=FALSE)

