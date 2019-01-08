library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

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


argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))


names(argsL) <- argsDF$V1


## Arg1 default
if(is.null(argsL$arg1)) {

	print ("Path to the Phenotype data,.csv file is missing")  
	quit()
}
 
## Arg2 default
if(is.null(argsL$arg2)) {
	print("Path to the data Directory is missing")
	quit()
}
 
## Arg3 default
if(is.null(argsL$arg3)) {
	print("Sample pattern is missing")
	quit()
}

if(is.null(argsL$arg4)) {
	print('covariate is missing')
	quit()
}

if(is.null(argsL$arg5)) {
	print (" adjustvars is missing")
	quit()
}

if(is.null(argsL$arg6)) {
	print ("Outfile for transcript is missing")
	quit()
}

if(is.null(argsL$arg7)) {
	print ("Outfile for gene is missing")
	quit()
}



pheno = read.csv(argsL$arg1)


bg = ballgown(dataDir=argsL$arg2, samplePattern = argsL$arg3, pData=pheno)
print ('ok')

bg_filt  = subset(bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

results_trans = stattest(bg_filt, feature="transcript", covariate=argsL$arg4, adjustvars=c(argsL$arg5), getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate=argsL$arg4, adjustvars=c(argsL$arg5), getFC=TRUE, meas="FPKM")

results_transcript = data.frame(geneNames=ballgown::geneNames(bg_filt),geneIDs=ballgown::geneIDs(bg_filt),results_trans)

results_transcript = arrange(results_transcript,pval)
results_genes = arrange(results_genes,pval)

write.csv(results_transcript,argsL$arg6, row.names=FALSE)
write.csv(results_genes, argsL$arg7, row.names=FALSE)

subset(results_genes,results_genes$qval<0.05)
subset(results_transcript,results_transcript$qval<0.05)

