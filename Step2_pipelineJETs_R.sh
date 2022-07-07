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
RprojectDir= _introducePATH_ #path to R
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



#-------------------------------------------------------
##-- STEP 2: Identification and classification of junctions, selection of JETs
#-------------------------------------------------------


################################################################################
##### 
##### Software path
RprojectDir= _introducePATH_ #path to R
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


	echo -e "Starting R -----" >> ${logFile}
	
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

		#specificaly calling STAR output files
		samFile="${prefix}_Chimeric.out.sam"
		bamFile="${prefix}_Aligned.sortedByCoord.out.bam"
		bamChimericFile="${prefix}_Chimeric.out.bam"
		bamChimericSortFile="${prefix}_Chimeric.out.sort.bam"
		chimericFile="${prefix}_Chimeric.out.junction"
		junctionFile="${prefix}_SJ.out.tab"

		logFinalOut="${prefix}_Log.final.out"
		libsize="$(grep "Uniquely mapped reads number" ${logFinalOut} |  sed -r 's/[\\t]+//g' | cut -d '|' -f2)" 

	

		echo -e "\e[1m${name}\t Calling and naming output files:\e[0m" >> ${logFile}
		echo -e "\e[1m${name}\t\t ${prefix}\e[0m" >> ${logFile}
		echo -e "\e[1m${name}\t\t ${libsize}\e[0m" >> ${logFile}


		size=11
	
			
			############################################################################
			##### Naming sample-specific output files
	
		fastaFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.fasta"
		idsFile="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.ids.txt"
		netmhcpanFile4="${prefix}_Fusions.annotatedchimJunc2e7.size${size}.netmhcpan4.0.txt"



			############################################################################
			##### R analysis
		
	
		cmd="${RscriptDir}/JET_analysis_filtered.R  \
			--chimeric ${chimericFile} \
			--junction ${junctionFile} \
			--genome ${genome} \
			--size ${size} \
			--libsize ${libsize} \
			--prefix ${prefix} \
			--verbose"
		

	      echo ${cmd} 
		

			date=`date +"%Y%m%d_%Hh%Mm%Ss"`
			echo -e "${name}\t$R_fusion_rna\t${rcmd}\t${date}" >> ${logFile}


done < ${infoFile}

echo -e "END" >> ${logFile}	