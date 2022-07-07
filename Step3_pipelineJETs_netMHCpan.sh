#!/bin/bash

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Ares Rocanin-Arjo


#-------------------------------------------------------
##-- Configuration
#-------------------------------------------------------


## Cluster Configuration
## =================================================================================================


#If necessary set up PATH environment (for instance: PATH=$PATH:/data/user/tools/CentOS/netMHCpan-3.0/) or :

##### Software path

netMHCpan4Dir= _introducePATH_ #path to NetMHCpan4

################################################################################
##### Set variables
day=`date +"%Y%m%d"`
time=`date +"%Hh%Mm%Ss"`
date="${day}_${time}"
readLength=100
organism="Mouse"
genome="mm10"
database="ensembl"


################################################################################
##### Configurate and set the following paths 

dataDir= _introducePATH_ #path to the input data directory
outputsDir= _introducePATH_ #path to the output results directory

logDir= _introducePATH_ #path to the logs directory

#path to the samplesInfo or metadata to perform a loop of submitted jobs
infoDir= _introducePATH_

#path to the star indexes. once done for always use the same
starIndexesDir= _introducePATH_

#path to the reference mm10 DONE (genetic or repeat masker)
metadataDir= _introducePATH_

#path to the tmp output and error files -- set up more if needed
ErrorDir= _introducePATH_


#if needed create those directories
mkdir -p ${logDir}

# path to the directory with the JET_identification_classification_Rscript.
RscriptDir= _introducePATH_



################################################################################
##### Calling DATA / creating files

	infoFile= _introduce PATH/filename_ #example "${infoDir}/Name_of_your_MetadataFile.txt" 
	

	logFile="${logDir}/R_netMHCpan_run_${date}.log"
	touch ${logFile}


################################################################################
##### Setting a loop when aligning several samples - optional

	while read rnaSample name day
	do
		echo -e "\e[1m${rnaSample}\t${name}\e[0m" >> ${logFile}

		
		############################################################################
		#####  Sample - specific Directories path -- if needed
		
		outputSampleDir="${outputsDir}/${name}_${day}"
		mkdir -p ${outputSampleDir} #making the directory if it do not exist

		echo -e "\e[1m${name}\t Creating the outputFiles and setting tmpDIR}\e[0m" >> ${logFile}

		

		############################################################################
		##### Calling sample-specific data files
	

		prefix="${outputSampleDir}/${name}" #setting up a prefix for each sample
	

		echo -e "\e[1m${name}\t Calling and naming output files:\e[0m" >> ${logFile}
		echo -e "\e[1m${name}\t\t ${prefix}\e[0m" >> ${logFile}
		echo -e "\e[1m${name}\t\t ${libsize}\e[0m" >> ${logFile}
	
			
			############################################################################
			##### Naming sample-specific output files
	
		fastaFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.fasta"
		idsFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.ids.txt"
		netmhcpanFile4="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.netmhcpan4.0.txt"

#-------------------------------------------------------
##-- STEP 3: Binding prediction to MHC class I
#-------------------------------------------------------
   
	echo -e "Starting size loop" >> ${logFile} #for each JET sequence we performed a peptide prediction of 8-11 aa

    	for size in `seq 8 11`
	 	  do

			############################################################################
			##### netMHCpan4.0

			"${netMHCpanDir4}/netMHCpan -a H-2-Kb,H-2-Db -f ${fastaFile} -l ${size} -xls -xlsfile ${netmhcpanFile4} -inptype 0 > ${netmhcpanFile4}.tmp"
			

			date=`date +"%Y%m%d_%Hh%Mm%Ss"`
			echo -e "${name}\tNetmhcpan\t${nmp4}\t${date}" >> ${logFile}

	done


done < ${infoFile}

echo -e "END" >> ${logFile}	