#library(ballgown)
#library(RSkittleBrewer)
#library(genefilter)
#library(dplyr)
#library(devtools)

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}
 
## Help section
if("--help" %in% args) {
  cat("
      The R Script
 
      Arguments:

      --arg1=someValue   - Path to the Phenotype data,.csv file 
      --arg2=someValue   - Path to the data Directory
      --arg3=someValue   - sample pattern 
      --arg4=someValue   - covariate 
      --arg5=someValue   - adjustvars
      --arg6=someValue   - Output .csv for transcript 
      --arg7=someValue   - Output .csv for Gene

      --help
 
      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")
 
  q(save="no")
}
 
## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

#print(parseArgs)

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))

#print (argsDF)
argsL <- as.list(as.character(argsDF$V2))


names(argsL) <- argsDF$V1


## Arg1 default
#if(is.null(args$arg1)) {
  ## do something
#}
 
## Arg2 default
#if(is.null(args$arg2)) {
  ## do something
#}
 
## Arg3 default
#if(is.null(args$arg3)) {
  ## do something
#}



"/home/joshij/Desktop/Sample_II/Rat_organs_phenotype_Data.csv"

'/home/joshij/Desktop/Sample_II/SE_result/Result_files/Merged_STRINGTIE_Out_files'

'SRR'



pheno = read.csv(argsL$arg1)

print (pheno)
bg = ballgown(dataDir=argsL$arg2, samplePattern = argsL$arg2, pData=pheno)
bg_filt  = subset(bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

results_trans = stattest(bg_filt, feature="transcript", covariate=argsL$arg4, adjustvars=c(argsL$arg5), getFC=TRUE, meas="FPKM")

results_genes = stattest(bg_filt, feature="gene", covariate=argsL$arg4, adjustvars=c(argsL$arg5), getFC=TRUE, meas="FPKM")

results_transcript = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),results_trans)
results_transcript = arrange(results_transcript,pval)
results_genes = arrange(results_genes,pval)

write.csv(results_transcript,'/home/joshij/Desktop/Sample_II/SE_result/Expression/chrX_transcript_results.csv', row.names=FALSE)
write.csv(results_genes, '/home/joshij/Desktop/Sample_II/SE_result/Expression/chrX_genes_results.csv', row.names=FALSE)

subset(results_genes,results_genes$qval<0.05)
ubset(results_transcript,results_transcript$qval<0.05)

