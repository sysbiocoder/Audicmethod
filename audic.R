#This is differential expression analysis for RNA sequencing data with no replicates using Audic and Claverie method 1997



# Input file
#The input file should comtain Genenames followed by counts and length columns separated by <TAB> space. 
#Note length column should be at the end.

# To run from command line 
# use  R CMD BATCH  '--args ctrl_column_no tmt_column_no input_file outfolder_name Length_column_no audic.R audic.Rout


require(stringr)

# Total number of reads
rm(list = ls())
# Read feature count 
	args = commandArgs() 
	ix = grep("--args", args, fixed = TRUE)
 	if (ix < length(args))
	{
    		for (i in (ix+1):length(args))
		 {
    	 		print(paste("argument #",i-ix,"-->",args[i], sep=" "), sep = "\n")  
		}
		#Parameters passed
		ctl<-args[8]
		ctl
		tmt <- args[9]
		tmt
		infile <- args[10]
		infile
		outfolder <- args[11]
		outfolder

		lenth <- args[12]
	} 
library(edgeR)
infile
ctl
tmt
data <- read.table(infile, header=TRUE,sep="\t")
outfolder
ctcol<-strtoi(ctl)
tmcol <- strtoi(tmt)
lencol <- strtoi(lenth)	
f1<-names(data)[ctcol]
f2<-names(data)[tmcol]
len <- names(data)[lencol]
ofile <- paste(outfolder,str_trim(f1),"_",str_trim(f2),sep="")
newdata <- data[,c(f1,f2,len)]
row.names(newdata) <- data$Geneid
N1 <- sum(newdata[[f1]])
N2 <- sum(newdata[[f2]])

x <- newdata[[f1]]
y <- newdata[[f2]]
newdata$csh <- newdata[[f1]]
newdata$tmt <- newdata[[f2]]

counts <- data.frame(newdata$csh,newdata$tmt)
dge <- DGEList(counts=counts,genes=data.frame(Length=newdata$Length))
dge <- calcNormFactors(dge)
RPKM <- rpkm(dge)
head(RPKM)

#------------------------------------Precalculations--------------------------------------------------------------------
d <- x+y
a <- factorial(x)*factorial(y)
b <- 1+(N2/N1)
c <- (x+y+1)
A <- (N2/N1)^y
C <- a*b^c
B <- factorial(x+y)


#----------------------------------------Audic----------------------------------------------------------------------------------

Pvalue <- exp(lchoose(x+y,x)-(x+y+1)*log(2))
FDR <-  p.adjust(Pvalue,method="bonferroni")
#---------------------FC----calculation-------------------------------------------------------------------------------------------

log2FC <- log2(RPKM[,1]/RPKM[,2])
#----------------------------------------Master file----------------------------------------------------------------------------------
acr <-data.frame(data$Geneid,log2FC,Pvalue,FDR)
file1 <- paste(ofile,"Master.txt",sep="")
write.table(acr,file1,sep="\t",row.names=F,quote=F)
