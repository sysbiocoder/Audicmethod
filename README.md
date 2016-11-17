# Audicmethod
RNAseq analysis without replicate
This is differential expression analysis for RNA sequencing data with no replicates using Audic and Claverie method 1997

To run from command line 

# Input file
The input file should be in the format containing Genenames, counts and length columns separated by <TAB> space

# usage 

R CMD BATCH  '--args ctrl_column_no tmt_column_no input_file outfolder_name Length_column_no audic.R audic.Rout
