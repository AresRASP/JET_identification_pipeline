#!/bin/bash

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Ares Rocanin-Arjo

#-------------------------------------------------------
##-- Configuration
#-------------------------------------------------------

#If necessary set up PATH envirovment or :

##### Software path
samtoolsBinDir= _introducePATH_ #path to samtools
starBinDir= _introducePATH_ #path to STAR-2.5.3a/bin/Linux_x86_64 


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


dataDir= _introducePATH_ #path to the input data dir
outputsDir= _introducePATH_ #path to the output results dir

logDir= _introducePATH_ #path to the logs dir

#path to the samplesinfo or metadata to perform a loop of submited jobs
infoDir= _introducePATH_

#path to the star indexes. once done for always use the same
starIndexesDir= _introducePATH_

#path to the refference mm10 DONE (genetic or repeatmasker)
metadataDir= _introducePATH_

#path to the tmp output and error files -- set up more if needed
ErrorDir= _introducePATH_


#if needed create those directories
mkdir -p ${logDir}



################################################################################
##### Calling DATA / creating files

	infoFile= _introduce PATH/filename_ #example "${infoDir}/Name_of_your_MetadataFile.txt" 
	
	#path to your refference fasta File
	fastaFile= _PATH/filename_FASTA_
	
	#path to your refference gtf File (genetic or repeatmasker)
	gtfGeneFile=  _PATH/filename_gtf_  # for instance "${metadataDir}/Mus_musculus.GRCm38.91.modified_oct18.gtf"  
	

	#Setting a script log - optional
	logFile="${logDir}/STARrun_${date}.log" #create a log File
	touch ${logFile}


#-------------------------------------------------------
##-- STEP 0: STAR pre-indexing
#-------------------------------------------------------

################################################################################
##### Pre-processing refference 
	threads=8
	cmd="${starBinDir}/STAR \
		--runThreadN ${threads} \
		--runMode genomeGenerate \
		--genomeDir ${starIndexesDir} \
		--genomeFastaFiles ${fastaFile} \
		--sjdbGTFfile ${gtfGeneFile} \
		--sjdbOverhang ${readLength}"
	echo ${cmd}
	#eval "time ${cmd}"


#-------------------------------------------------------
##-- STEP 1: STAR alignment
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

	
	while read rnaSample name #example: rnaSample="D804T348"; name="Thymus_4"
	do
		echo -e "\e[1m${rnaSample}\t${name}\e[0m" >> ${logFile}



		############################################################################
		#####  Sample - specific Directories path -- if needed

		outputSampleDir="${outputsDir}/${name}_${day}" 
		mkdir -p ${outputSampleDir}

		#tmp dir to save tmp STAR files 
		tmpDir=
		
		echo -e "\e[1m${rnaSample}\t Creating the outputFiles and setting tmpDIR}\e[0m" >> ${logFile}



		############################################################################
		##### Calling sample-specific data and naming output files
		
		# path to sampleFastq
	   	fastqR1Filegz= #example: "${dataDir}/${rnaSample}.R1.fastq.gz"
		fastqR2Filegz= #example: "${dataDir}/${rnaSample}.R2.fastq.gz"

		echo -e "\e[1m${rnaSample}\t Reading fasta files\t${fastqR1File} i ${fastqR2File}\e[0m" >> ${logFile}
		

		prefix="${outputSampleDir}/${name}" #setting up a prefix for each sample

		#specificaly naming output files
		samFile="${prefix}_Chimeric.out.sam"
		bamFile="${prefix}_Aligned.sortedByCoord.out.bam"
		bamChimericFile="${prefix}_Chimeric.out.bam"
		bamChimericSortFile="${prefix}_Chimeric.out.sort.bam"
		chimericFile="${prefix}_Chimeric.out.junction"
		junctionFile="${prefix}_SJ.out.tab"


		echo -e "\e[1m${rnaSample}\t Naming output files:\e[0m" >> ${logFile}
		echo -e "\e[1m${rnaSample}\t\t ${prefix}\e[0m" >> ${logFile}



		############################################################################
		##### Test if files exist
		
		flag=false

      	# OPTIONAL: decompress fastq files
		if [ ! -f ${fastqR1File} -a -f ${fastqR1Filegz} ] ; then gunzip -c $fastqR1Filegz > $fastqR1File; fi
		if [ ! -f ${fastqR2File} -a -f ${fastqR2Filegz} ] ; then gunzip -c $fastqR2Filegz > $fastqR2File; fi

		#Cheking the sample fastq files exist
		if [ ! -f ${fastqR1File} -a ! -f ${fastqR1Filegz} ] ; then echo "${fastqR1File} not found !" 1>&2 ; fastqR1File="NA" ; flag=true ; fi
		if [ ! -f ${fastqR2File} -a ! -f ${fastqR2Filegz} ] ; then echo "${fastqR2File} not found !" 1>&2 ; fastqR2File="NA" ; flag=true ; fi
	
		
		if ${flag} ; then continue; fi
		echo -e "${name}\t${fastqR1File}\t${fastqR2File}" >> ${logFile}
		


		############################################################################
		##### STAR alignment configuration
	
		threads=8
		cmd="${starBinDir}/STAR \
				--quantMode GeneCounts \
				--twopassMode Basic \
				--runThreadN ${threads} \
				--genomeDir ${starIndexesDir} \
				--sjdbGTFfile ${gtfGeneFile} \
				--sjdbOverhang 100 \
				--readFilesIn ${fastqR1File} ${fastqR2File} \
				--outFileNamePrefix ${prefix}_ \
				--outTmpDir ${tmpDir}/STAR_${name}_${date} \
				--outReadsUnmapped Fastx \
				--outSAMtype BAM SortedByCoordinate \
				--bamRemoveDuplicatesType UniqueIdentical \
				--outFilterMismatchNoverLmax 0.04 \
				--outMultimapperOrder Random \
				--outFilterMultimapNmax 1000 \
				--winAnchorMultimapNmax 1000 \
				--chimOutType WithinBAM \
				--chimSegmentMin 10 \
				--chimJunctionOverhangMin 10 ; \
			${samtoolsBinDir}/samtools view -@ ${threads} -b ${samFile} > ${bamChimericFile} ; \
			${samtoolsBinDir}/samtools sort -@ ${threads} -o ${bamChimericSortFile} -O bam ${bamChimericFile} ; \
			${samtoolsBinDir}/samtools index ${bamChimericSortFile} ; \
			${samtoolsBinDir}/samtools index ${bamFile}"
		
		echo ${cmd}

		date=`date +"%Y%m%d_%Hh%Mm%Ss"`
		echo -e "${name}\tStar_rna:\t${starcmd}\t${date}" >> ${logFile}


done < ${infoFile}
