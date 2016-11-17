#RNAseq analysis without replicate using Audicmethod
This program can be used to identify the differential gene expression analysis of raw counts from RNA sequencing data without replicates using Audic and Claverie method 1997

# Input file
The input file should be in the format containing Genenames, counts and length columns separated by <TAB> space

 
# usage 
To run from command line
R CMD BATCH  '--args ctrl_column_no tmt_column_no input_file outfolder_name Length_column_no' audic.R audic.Rout
