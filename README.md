# README for "JET identification pipeline" repository

-------------------------------------------------------
###    Created by Alexandre Houy and Christel Goudot
###    Modified and adapted by Ares Rocanin-Arjo
###    				07/2022
-------------------------------------------------------

This pipeline was used in the manuscript: Epigenetically-controlled tumor antigens derived from splice junctions between exons and transposable elements (Burbage, M. and Rocanin-Arjo, A.) to identify, quantify and select non-canonical splice junctions between exons and transposable elements (or JETs).
It consists in 3 bash (.sh) scripts and one Rscript.  
	1- Step1_pipelineJETs_STAR.sh  
	2- Step2_pipelineJETs_R.sh  
	3- JET_analysis_filtered.R  
	4- Step3_pipelineJETs_netMHCpan.sh  

The R script, is called by one of the sh files(the Step2). 
In all three of them there are paths and files that must be defined to adapt the script to your computer environment and references files. This are indicated as `_introducePATH_` or `_introduce PATH/filename_` (also `_PATH/filename_gtf_` or `_PATH/filename_FASTA_`) in the scipts.


This particular script requires the following softwares and versions:
- STAR 2.5.3
- samtools-1.3
- Rproject software 3.2.3 with packages biovizBase1.18.0, GenomicAligments1.6.3, GenomicFeatures_1.22.13, data.table_1.12.8, scales_1.1.0, ggplot2_3.2.1.

It will also require a repeat masker annotation file from USCS. The indication are:
#### INDICATIONS TO DOWNLOAD ANNOTATION: 
Files downloaded from UCSC table --- done once  

- Go to : Tools - Table browser  
    * clade               : Mammal  
    * genome              : <organism>  
    * assembly            : <genome>  
    * group               : Variation and Repeats  
    * track               : RepeatMasker  
    * table               : rmsk  
    * region              : genome  
    * output format       : select fields from primary and related table  
    * output file         : <genome>.repeatmasker_reformat.txt  
    * file type returned : plain text  
  -> Click get output  
	
 - Select  
		* genoName  
		* genostart  
		* genoend  
		* strand  
		* repName  
		* repClass  
		* repfamily  
 -> Click get output  

Also the path and filename must be indicated in the R script in the 118th line.  
  
### Citation
Please use the following information to cite:

### Contact the Author
Ares Rocanin-Arjo: maria-ares.rocanin-arjo@curie.fr, Marianne Burbage: marianne.burbage@curie.fr and Christel Goudot: christel.goudot@curie.fr
