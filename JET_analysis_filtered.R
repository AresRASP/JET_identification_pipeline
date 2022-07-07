#!/data/users/mrocanin/R/R-3.2.3/bin/Rscript

#### created by Alexandre Houy and Christel Goudot
#### modified and adapted by Ares Rocanin-Arjo

rm(list = ls(all.names = T))
invisible(gc())
graphics.off()

tall = Sys.time()


#!!!
r.path    = _introducePATH_  # Local directory of R packages
R         = paste(R.Version()$major, strsplit(x = as.character(R.Version()$minor), split = "\\.")[[1]][1], sep = ".")
#!!!
lib.local = paste( _instroducePATH_ ) #path to R library if specific

#dir.create(lib.local, showWarnings = F, recursive = T)
.libPaths(new = lib.local)


options(useHTTPS = F, BioC_mirror="http://bioconductor.org")




##### Getops #####
#--------------------------------------------------------------------------------

suppressMessages(require(getopt))
spec = matrix(c("chimeric", "c", 1, "character",
                "junction", "j", 1, "character",
                "genome",   "g", 1, "character",
                "prefix",   "p", 1, "character",
                "size",     "s", 1, "integer",
                "libsize",	"l", 1, "integer",
                "help",     "h", 0, "logical",
                "verbose",  "v", 0, "logical"),
              byrow = TRUE,
              ncol = 4);
opt  = getopt(spec)

if(! is.null(opt$help)) stop(getopt(spec, usage = T))

verbose        = ifelse(is.null(opt$verbose), F, T)           # Verbose mode
size           = ifelse(is.null(opt$size), 9, opt$size)       # Peptide size
chimeric.file  = opt$chimeric                                 # Chimeric.out.junction file from STAR
junctions.file = opt$junction                                 # SJ.out.tab file from STAR
prefix         = opt$prefix                                   # Prefix of output files & plots
genome         = opt$genome                                   # Genome version
libsize   	   = ifelse(is.null(opt$libsize), NA, opt$libsize) # libsize of the sample to normalize counts


## cheking if the data exists


flag = F
if(is.null(chimeric.file) & is.null(junctions.file)){             cat("/!\\ Error : At least one parameter (--chimeric or --junction) must be set /!\\\n", file = stderr()); flag = T }
if(! is.null(chimeric.file)){ if(! file.exists(chimeric.file)){   cat("/!\\ Error : Argument --chimeric", chimeric.file, "not found /!\\\n", file = stderr());               flag = T } }
if(! is.null(junctions.file)){ if(! file.exists(junctions.file)){ cat("/!\\ Error : Argument --junction", junctions.file, "not found /!\\\n", file = stderr());              flag = T } }
if(! genome %in% c("hg19", "mm9", "mm10")){                       cat("/!\\ Error : Argument --genome", genome, "must be 'hg19', 'mm9' or 'mm10' /!\\\n", file = stderr());  flag = T }
if(! is.integer(libsize)){                                        cat("/!\\ Error : Argument --libsize must be an integer /!\\\n", file = stderr());                         flag = T }                                                   
if(flag) stop(getopt(spec, usage = T))




##### Libraries #####
#--------------------------------------------------------------------------------


if(verbose) cat("Loading packages")
t = Sys.time()
ipak = function(list.packages){
	list.new.packages =  list.packages[!(list.packages %in% installed.packages()[, "Package"])]
	if (length(list.new.packages)){
		suppressMessages(source("http://bioconductor.org/biocLite.R"))
		biocLite(list.new.packages, dependencies = TRUE)
	}
	sapply(list.packages, require, character.only = T, quietly = T)
}



list.packages = c("GenomicFeatures",
                  "GenomicAlignments",
                  "data.table",
                  "scales",
                  "ggbio",
                  "ggplot2",
                  "biovizBase")
ipak.res = ipak(list.packages = list.packages)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")







##### Calling refference annotations #####
#--------------------------------------------------------------------------------

if(genome %in% c("hg19")){
	organism    = "Human"
	chromosomes = paste0("chr", c(1:22, "X"))
	bs          = paste("BSgenome", "Hsapiens", "UCSC", genome, sep = ".")
	txdb        = paste("TxDb", "Hsapiens", "UCSC", genome, "ensGene", sep = ".")
} else if(genome %in% c("mm9", "mm10")){
	organism    = "Mouse"
	chromosomes = paste0("chr", c(1:19, "X"))
	bs          = paste("BSgenome", "Mmusculus", "UCSC", genome, sep = ".")
	txdb        = paste("TxDb", "Mmusculus", "UCSC", genome, "ensGene", sep = ".")
}

#!!!
repeats.file = _introduce PATH/name_ #path to the reeatmasker obtained at UCSC as indicated in the README

if(verbose) cat("Loading genome data")

t = Sys.time()
data(ideo)
genome.ideo         = suppressMessages(ideo[[genome]])
genome.ideo         = keepSeqlevels(x = genome.ideo, value = chromosomes ) #ADDING HERE THE OPTION PRUNING.MODE coarse! error option by default does no pruning and results with error in the "pipeline" #pruning.mode="coarse"
genome(genome.ideo) = genome
	
ipak.gen            = ipak(c(txdb, bs))
genome.bs           = eval(parse(text = bs))
genome.txdb         = eval(parse(text = txdb))
genome.txdb         = keepSeqlevels(x = genome.txdb, value = chromosomes)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")




##### Setting Functions #####
#--------------------------------------------------------------------------------

get.overlap = function(query, subject, type = "any", name = NULL){
	hits                 = as.data.frame(findOverlaps(query = query, subject = subject, type = type, select = "all"))
	hits                 = as.data.table(data.frame(queryHits = hits$queryHits, as.data.frame(subject[hits$subjectHits], row.names = NULL)))
	hits                 = hits[, lapply(.SD, function(x) paste(sort(as.character(unlist(x))), collapse = ";")), by = queryHits]
	if(! is.null(name)) colnames(hits) = paste(colnames(hits), name, sep = ".")
	df                   = as.data.frame(matrix(data = NA, nrow = length(query), ncol = ncol(hits) - 1, dimnames = list(c(NULL), c(colnames(hits)[-1]))))
	hits=data.frame(hits)
	df[hits$queryHits, ] = hits[, -1]
	return(df)
}

build.fusion.exon = function(type = c("donor.end", "acceptor.start"), s.gr, tx.id){
	# Adapted from 'chimera' package (function '.BuildFusion')
	tx.list = list("tx_id" = tx.id)
	e.gr    = exons(genome.txdb, vals = tx.list, columns = c("tx_id", "tx_name", "exon_rank"))
	f.pos   = findOverlaps(query = s.gr, subject = e.gr, type = "any", select = "first", ignore.strand = T)
	e.tx    = as.list(elementMetadata(e.gr)$tx_id)
	e.rank  = elementMetadata(e.gr)$exon_rank
	e.idx   = sapply(e.tx, function(x, y){ grep(y, x) }, y = tx.id)
	rank.e  = NULL
	for(i in 1:length(e.idx)){ rank.e[i] = e.rank[[i]][e.idx[i]] }
	if(type == "donor.end"){
		if(!is.na(f.pos)){
			if(unique(as.character(strand(e.gr))) == "-"){
				start(e.gr[f.pos]) = start(s.gr)
				e.gr  = e.gr[f.pos:length(e.gr)]
				e.seq = getSeq(genome.bs, e.gr)
				e.rnk = seq(length(e.seq), 1)
			}
			else{
				end(e.gr[f.pos]) = end(s.gr)
				e.gr  = e.gr[1:f.pos]
				e.seq = getSeq(genome.bs, e.gr)
				e.rnk = seq(1, length(e.seq))
			}
			donor.seq  = NULL
			for(i in e.rnk){ donor.seq = c(donor.seq, as.character(e.seq[i])) }
			donor.seq = paste(donor.seq, collapse = "")
			donor.seq = DNAStringSet(donor.seq)
			return(donor.seq)
		}
		else{
			return(NULL)
		}
	}
	else if(type == "acceptor.start"){
		if(!is.na(f.pos)){
			if(unique(as.character(strand(e.gr))) == "-"){
				end(e.gr[f.pos]) = end(s.gr)
				e.gr  = e.gr[1:f.pos]
				e.seq = getSeq(genome.bs, e.gr)
				e.rnk = seq(length(e.seq), 1)
			}
			else{
				start(e.gr[f.pos]) = start(s.gr)
				e.gr  = e.gr[f.pos:length(e.gr)]
				e.seq = getSeq(genome.bs, e.gr)
				e.rnk = seq(1, length(e.seq))
			}
			acceptor.seq = NULL
			for(i in e.rnk){ acceptor.seq = c(acceptor.seq, as.character(e.seq[i])) }             
			acceptor.seq = paste(acceptor.seq, collapse = "")
			acceptor.seq = DNAStringSet(acceptor.seq)
			return(acceptor.seq)
		}
		else{
			return(NULL)
		}
	}
}

build.fusion.repeat = function(type = c("donor.end", "acceptor.start"), s.gr, r.gr){
	if(type == "donor.end"){
		if(unique(as.character(strand(r.gr))) == "-"){ start(r.gr) = start(s.gr) }
		else{ end(r.gr) = end(s.gr) }
		donor.seq = getSeq(genome.bs, r.gr)
		return(donor.seq)
	}
	else if(type == "acceptor.start"){
		if(unique(as.character(strand(r.gr))) == "-"){ end(r.gr) = end(s.gr) }
		else{ start(r.gr) = start(s.gr) }
		acceptor.seq = getSeq(genome.bs, r.gr)
		return(acceptor.seq)
	}
}

get.fusion.sequence = function(donor, acceptor){
	donor.seq = NULL
	if(! is.na(donor$tx_id.exon) & is.na(acceptor$tx_id.exon)){
		tx_ids   = unlist(strsplit(x = donor$tx_id.exon, split = ";"))
		tx_names = unlist(strsplit(x = donor$tx_name.exon, split = ";"))
		for(i in 1:length(tx_ids)){
			donor.seq.tmp        = build.fusion.exon(type = "donor.end", s.gr = donor, tx.id = tx_ids[i])
			names(donor.seq.tmp) = paste(seqnames(donor), start(donor), strand(donor), tx_names[i], sep = ":")
			donor.seq            = c(donor.seq, donor.seq.tmp)
		}
	}else if(! is.na(donor$superfamily.repeat)){
		donor.seqnames    = unlist(strsplit(x = donor$seqnames.repeat, split = ";"))
		donor.start       = unlist(strsplit(x = donor$start.repeat, split = ";"))
		donor.end         = unlist(strsplit(x = donor$end.repeat, split = ";"))
		donor.strand      = unlist(strsplit(x = donor$strand.repeat, split = ";"))
		donor.superfamily = unlist(strsplit(x = donor$superfamily.repeat, split = ";"))
		donor.family      = unlist(strsplit(x = donor$family.repeat, split = ";"))
		donor.subfamily   = unlist(strsplit(x = donor$subfamily.repeat, split = ";"))

		for(i in 1:length(donor.seqnames)){
			donor.tmp            = with(donor, GRanges(seqnames    = donor.seqnames[i],
			                                           ranges      = IRanges(start = as.numeric(donor.start[i]), end = as.numeric(donor.end[i])),
			                                           strand      = donor.strand[i],
			                                           superfamily = donor.superfamily[i],
			                                           family      = donor.family[i],
			                                           subfamily   = donor.subfamily[i]))
			donor.seq.tmp        = build.fusion.repeat("donor.end", s.gr = donor, r.gr = donor.tmp)
			names(donor.seq.tmp) = paste(seqnames(donor), start(donor), strand(donor), donor.tmp$subfamily[i], sep = ":")
			donor.seq            = c(donor.seq, donor.seq.tmp)
		}
	}

	acceptor.seq = NULL
	if(! is.na(acceptor$tx_id.exon) & is.na(donor$tx_id.exon)){
		tx_ids   = unlist(strsplit(x = acceptor$tx_id.exon, split = ";"))
		tx_names = unlist(strsplit(x = acceptor$tx_name.exon, split = ";"))
		for(i in 1:length(tx_ids)){
			acceptor.seq.tmp        = build.fusion.exon(type = "acceptor.start", s.gr = acceptor, tx.id = tx_ids[i])
			names(acceptor.seq.tmp) = paste(seqnames(acceptor), start(acceptor), strand(acceptor), tx_names[i], sep = ":")
			acceptor.seq            = c(acceptor.seq, acceptor.seq.tmp)
		}
	}else if(! is.na(acceptor$superfamily.repeat)){
		acceptor.seqnames    = unlist(strsplit(x = acceptor$seqnames.repeat, split = ";"))
		acceptor.start       = unlist(strsplit(x = acceptor$start.repeat, split = ";"))
		acceptor.end         = unlist(strsplit(x = acceptor$end.repeat, split = ";"))
		acceptor.strand      = unlist(strsplit(x = acceptor$strand.repeat, split = ";"))
		acceptor.superfamily = unlist(strsplit(x = acceptor$superfamily.repeat, split = ";"))
		acceptor.family      = unlist(strsplit(x = acceptor$family.repeat, split = ";"))
		acceptor.subfamily   = unlist(strsplit(x = acceptor$subfamily.repeat, split = ";"))

		for(i in 1:length(acceptor.seqnames)){
			acceptor.tmp             = with(acceptor, GRanges(seqnames    = acceptor.seqnames[i],
			                                                 ranges      = IRanges(start = as.numeric(acceptor.start[i]), end = as.numeric(acceptor.end[i])),
			                                                 strand      = acceptor.strand[i],
			                                                 superfamily = acceptor.superfamily[i],
			                                                 family      = acceptor.family[i],
			                                                 subfamily   = acceptor.subfamily[i]))
			acceptor.seq.tmp        = build.fusion.repeat("acceptor.start", s.gr = acceptor, r.gr = acceptor.tmp)
			names(acceptor.seq.tmp) = paste(seqnames(acceptor), start(acceptor), strand(acceptor), acceptor.tmp$subfamily[i], sep = ":")
			acceptor.seq            = c(acceptor.seq, acceptor.seq.tmp)
		}
	}

	sequences = NULL
	for(i in 1:length(donor.seq)){
		for(j in 1:length(acceptor.seq)){
			sequence.tmp = list(list(Donor = donor.seq[i][[1]], Acceptor = acceptor.seq[j][[1]]))
			names(sequence.tmp) = paste(names(donor.seq[i][[1]]), names(acceptor.seq[j][[1]]), sep = ">")
			sequences = c(sequences, sequence.tmp)
		}
	}
	return(sequences)
}



get.peptides = function(sequences, ids, order = c("R1R2", "R2R1")){
	name    = names(sequences)

	for(j in 1:length(name)){
		donor.full      = sequences[[name[j]]]$Donor
		donor.substr    = substr(x = donor.full, start = nchar(donor.full) + 1 - (size * 3 - 1), stop = nchar(donor.full))
		acceptor.full   = sequences[[name[j]]]$Acceptor
		acceptor.substr = substr(x = acceptor.full, 1, size * 3 - 1)

		fusion.dna      = paste0(donor.substr, acceptor.substr)

		fusion.ORF1     = translate(DNAStringSet(substr(x = fusion.dna, start = 1, stop = nchar(fusion.dna))))
		fusion.ORF1     = gsub(pattern = "\\*.*", replacement = "", x = fusion.ORF1)
		name.ORF1       = paste(order, "ORF1", name[j], paste0("width=", nchar(fusion.ORF1)), sep = ";")
		if(fusion.ORF1 %in% names(ids)){ ids[[fusion.ORF1]] = c(ids[[fusion.ORF1]], name.ORF1) } else { ids[[fusion.ORF1]] = list(name.ORF1) }

		fusion.ORF2     = translate(DNAStringSet(substr(x = fusion.dna, start = 2, stop = nchar(fusion.dna))))
		fusion.ORF2     = gsub(pattern = "\\*.*", replacement = "", x = fusion.ORF2)
		name.ORF2       = paste(order, "ORF2", name[j], paste0("width=", nchar(fusion.ORF2)), sep = ";")
		if(fusion.ORF2 %in% names(ids)){ ids[[fusion.ORF2]] = c(ids[[fusion.ORF2]], name.ORF2) } else { ids[[fusion.ORF2]] = list(name.ORF2) }

		fusion.ORF3     = translate(DNAStringSet(substr(x = fusion.dna, start = 3, stop = nchar(fusion.dna))))
		fusion.ORF3     = gsub(pattern = "\\*.*", replacement = "", x = fusion.ORF3)
		name.ORF3       = paste(order, "ORF3", name[j], paste0("width=", nchar(fusion.ORF3)), sep = ";")
		if(fusion.ORF3 %in% names(ids)){ ids[[fusion.ORF3]] = c(ids[[fusion.ORF3]], name.ORF3) } else { ids[[fusion.ORF3]] = list(name.ORF3) }
		
	}

	return(ids)
}



 
##### Calling Data Files #####
#--------------------------------------------------------------------------------


if(verbose) cat("Read ", repeats.file, sep = "")

t = Sys.time()
repeats                = read.table(file             = repeats.file,
                                    sep              = "\t",
                                    stringsAsFactors = F,
                                    header           = F,
                                    col.names        = c("chr", "start", "end", "strand", "subfamily", "superfamily", "family"))
repeats$Id             = rownames(repeats)
repeats                = subset(repeats, chr %in% chromosomes)
repeats$start          = repeats$start + 1
repeats.gr             = makeGRangesFromDataFrame(repeats, keep.extra.columns = T)
seqlengths(repeats.gr) = seqlengths(genome.ideo)[names(seqlengths(repeats.gr))]
genome(repeats.gr)     = genome
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")


if(verbose) cat("Get exon positions")
t = Sys.time()
exons.gr = cds(genome.txdb, columns = c("gene_id", "tx_id", "tx_name", "exon_rank"))
exons.gr = keepSeqlevels(x = exons.gr, value = chromosomes)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")



if(verbose) cat("Get promoter positions")
t = Sys.time()
promoters.gr = promoters(genome.txdb, columns = c("gene_id", "tx_id", "tx_name", "exon_rank"), upstream = 2000, downstream = 2000)
promoters.gr = keepSeqlevels(x = promoters.gr, value = chromosomes)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

chimeric = NULL
if(! is.null(chimeric.file)){
	if(verbose) cat("Read", chimeric.file)
	t = Sys.time()
	chimeric          = read.table(file             = chimeric.file,
	                               sep              = "\t",
	                               stringsAsFactors = F,
	                               col.names        = c("chr.from", "pos.from", "str.from", "chr.to", "pos.to", "str.to",
	                                                    "type", "repeat.from", "repeat.to", "Name",
	                                                    "base.from", "CIGAR.from", "base.to", "CIGAR.to"))
	chimeric          = subset(chimeric, chr.from %in% chromosomes & chr.to %in% chromosomes)
	chimeric          = subset(chimeric, type %in% c(1, 2))
	chimeric$n        = 1
	chimeric          = aggregate(formula = n ~ chr.from + pos.from + str.from + repeat.from + chr.to + pos.to + str.to + repeat.to + type,
	                              data    = chimeric,
	                              FUN     = sum)
	

	########Expression filter
 	if(verbose) cat(format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "\tlibrary size filtering \n")
   	chimeric          = subset(chimeric, n/libsize >= 2*1e-7)
 
	


	chimeric$pos.from = chimeric$pos.from - 1
	chimeric$pos.to   = chimeric$pos.to + 1
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")
}



junctions = NULL
if(! is.null(junctions.file)){
	if(verbose) cat("Read", junctions.file)
	t = Sys.time()
	junctions          = read.table(file             = junctions.file,
	                               sep              = "\t",
	                               stringsAsFactors = F,
	                               col.names        = c("chr", "start", "end", "strand", "motif", "annotated", "unique", "multi", "overhang"))
	junctions          = subset(junctions, chr %in% chromosomes)
	junctions          = subset(junctions, motif %in% c(1, 2))
	junctions$strand   = with(junctions, ifelse(strand == 1, "+", ifelse(strand == 2, "-", "*")))
	junctions$pos.from = with(junctions, ifelse(strand == "+", start - 1, end + 1))
	junctions$pos.to   = with(junctions, ifelse(strand == "+", end + 1, start - 1))
	junctions          = with(junctions, data.frame(chr.from = chr, pos.from = pos.from, str.from = strand, repeat.from = 0,
	                                            chr.to = chr, pos.to = pos.to, str.to = strand, repeat.to = 0,
	                                            type = motif, n = unique + multi))
	

	#######Expression filter

 	 if(verbose) cat(format(x = Sys.time(), format = "%Y-%m-%d %H:%M:%S"), "\tlibrary size filtering \n")
   	 junctions          = subset(junctions, n/libsize >= 2*1e-7)

	

	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")
}



if(verbose) cat("Bind chimeric and junction tables")
t = Sys.time()
chimeric.all = rbind(data.frame(file = "Chimeric", chimeric),
                     data.frame(file = "Junction", junctions))
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")




	##### Donor #####
	#--------------------------------------------------------------------------------

	if(verbose) cat("Create GRanges (donor)")
	t = Sys.time()
	chimeric.from.gr             = with(chimeric.all, GRanges(seqnames = chr.from,
	                                                          ranges   = IRanges(start = pos.from, width = 1),
	                                                          strand   = str.from,
	                                                          file     = file,
	                                                          type     = type,
	                                                          n        = n))


	seqlengths(chimeric.from.gr) = seqlengths(genome.ideo)[names(seqlengths(chimeric.from.gr))]
	genome(chimeric.from.gr)     = genome
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")


	if(verbose) cat("Annotate repeat positions (donor)")
	t = Sys.time()
	repeats.from.tmp        = get.overlap(query = chimeric.from.gr, subject = repeats.gr, name = "repeat")
	mcols(chimeric.from.gr) = DataFrame(mcols(chimeric.from.gr), repeats.from.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")


	if(verbose) cat("Annotate gene positions (donor)")
	t = Sys.time()
	exons.from.tmp          = get.overlap(query = chimeric.from.gr, subject = exons.gr, name = "exon")
	mcols(chimeric.from.gr) = DataFrame(mcols(chimeric.from.gr), exons.from.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")


	if(verbose) cat("Annotate promoter positions (donor)")
	t = Sys.time()
	promoters.from.tmp      = get.overlap(query = chimeric.from.gr, subject = promoters.gr, name = "Promoter")
	mcols(chimeric.from.gr) = DataFrame(mcols(chimeric.from.gr), promoters.from.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")



	##### Acceptor #####
	#--------------------------------------------------------------------------------

	if(verbose) cat("Create GRanges (acceptor)")
	t = Sys.time()
	chimeric.to.gr             = with(chimeric.all, GRanges(seqnames = chr.to,
	                                                        ranges   = IRanges(start = pos.to, width = 1),
	                                                        strand   = str.to))
	seqlengths(chimeric.to.gr) = seqlengths(genome.ideo)[names(seqlengths(chimeric.to.gr))]
	genome(chimeric.to.gr)     = genome
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

	if(verbose) cat("Annotate repeat positions (acceptor)")
	t = Sys.time()
	repeats.to.tmp        = get.overlap(query = chimeric.to.gr, subject = repeats.gr, name = "repeat")
	mcols(chimeric.to.gr) = DataFrame(mcols(chimeric.to.gr), repeats.to.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

	if(verbose) cat("Annotate gene positions (acceptor)")
	t = Sys.time()
	exons.to.tmp          = get.overlap(query = chimeric.to.gr, subject = exons.gr, name = "exon")
	mcols(chimeric.to.gr) = DataFrame(mcols(chimeric.to.gr), exons.to.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

	if(verbose) cat("Annotate promoter positions (acceptor)")
	t = Sys.time()
	promoters.to.tmp      = get.overlap(query = chimeric.to.gr, subject = promoters.gr, name = "Promoter")
	mcols(chimeric.to.gr) = DataFrame(mcols(chimeric.to.gr), promoters.from.tmp)
	if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")



##### Combine and filter #####
#--------------------------------------------------------------------------------

if(verbose) cat("Combine and filter data")
t = Sys.time()
chimeric.combine.gr                         = chimeric.from.gr
values(chimeric.combine.gr)$acc             = chimeric.to.gr
chimeric.combine.gr                         = keepSeqlevels(chimeric.combine.gr, chromosomes)
seqlengths(chimeric.combine.gr)             = seqlengths(genome.ideo)[names(seqlengths(chimeric.combine.gr))]
chimeric.combine.gr$acc                     = keepSeqlevels(chimeric.combine.gr$acc, chromosomes)
seqlengths(chimeric.combine.gr$acc)         = seqlengths(genome.ideo)[names(seqlengths(chimeric.combine.gr$acc))]
filter.gene                                 = (! is.na(chimeric.to.gr$superfamily.repeat)   & ! is.na(chimeric.from.gr$seqnames.exon)     & is.na(chimeric.to.gr$seqnames.exon)) |
                                              (! is.na(chimeric.from.gr$superfamily.repeat) & is.na(chimeric.from.gr$seqnames.exon)       & ! is.na(chimeric.to.gr$seqnames.exon))
filter.promoter                             = (! is.na(chimeric.to.gr$superfamily.repeat)   & ! is.na(chimeric.from.gr$seqnames.Promoter) & is.na(chimeric.to.gr$seqnames.Promoter)) |
                                              (! is.na(chimeric.from.gr$superfamily.repeat) & is.na(chimeric.from.gr$seqnames.Promoter)   & ! is.na(chimeric.to.gr$seqnames.Promoter))
values(chimeric.combine.gr)$filter.gene     = filter.gene
values(chimeric.combine.gr)$filter.promoter = filter.promoter
values(chimeric.combine.gr)$filter          = filter.gene | filter.promoter
values(chimeric.combine.gr)$id              = 1:length(chimeric.combine.gr)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")





##### Table #####
#--------------------------------------------------------------------------------

table.file = paste0(prefix, "_Chimeric.out.annotatedJET.txt")
t = Sys.time()
if(verbose) cat("Write", table.file)
chimeric.df = as.data.frame(chimeric.combine.gr)
chimeric.df$Donor = with(chimeric.df, paste(ifelse(! is.na(seqnames.repeat), "R", "."),
                                            ifelse(! is.na(seqnames.exon), "E", "."),
                                            ifelse(! is.na(seqnames.Promoter), "P", "."), sep = ";"))
chimeric.df$Acceptor = with(chimeric.df, paste(ifelse(! is.na(acc.seqnames.repeat), "R", "."),
                                               ifelse(! is.na(acc.seqnames.exon), "E", "."),
                                               ifelse(! is.na(acc.seqnames.Promoter), "P", "."), sep = ";"))
chimeric.df = chimeric.df[, c("id", "file", "type", "n", "Donor", "Acceptor", "filter.gene", "filter.promoter", "filter",
                              "seqnames", "start", "end", "width", "strand",
                              "seqnames.repeat", "start.repeat", "end.repeat", "width.repeat", "strand.repeat", "subfamily.repeat", "superfamily.repeat", "family.repeat", "Id.repeat",
                              "seqnames.exon", "start.exon", "end.exon", "width.exon", "strand.exon", "gene_id.exon", "tx_id.exon", "tx_name.exon", "exon_rank.exon",
                              "seqnames.Promoter", "start.Promoter", "end.Promoter", "width.Promoter", "strand.Promoter", "gene_id.Promoter", "tx_id.Promoter", "tx_name.Promoter", "exon_rank.Promoter",
                              "acc.seqnames", "acc.start", "acc.end", "acc.width", "acc.strand",
                              "acc.seqnames.repeat", "acc.start.repeat", "acc.end.repeat", "acc.width.repeat", "acc.strand.repeat", "acc.subfamily.repeat", "acc.superfamily.repeat", "acc.family.repeat", "acc.Id.repeat",
                              "acc.seqnames.exon", "acc.start.exon", "acc.end.exon", "acc.width.exon", "acc.strand.exon", "acc.gene_id.exon", "acc.tx_id.exon", "acc.tx_name.exon", "acc.exon_rank.exon",
                              "acc.seqnames.Promoter", "acc.start.Promoter", "acc.end.Promoter", "acc.width.Promoter", "acc.strand.Promoter", "acc.gene_id.Promoter", "acc.tx_id.Promoter", "acc.tx_name.Promoter", "acc.exon_rank.Promoter")]
write.table(file = table.file, x = chimeric.df, sep = "\t", quote = F, row.names = F, col.names = T)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")





##### Extract sequences : JETs #####
#--------------------------------------------------------------------------------


if(verbose) cat("Get fusion sequences\n")
t = Sys.time()
chimeric.combine.gr.tmp = subset(chimeric.combine.gr, filter)
ids                     = NULL
sequences.all           = NULL
for(i in 1:length(chimeric.combine.gr.tmp)){
	# Read1 > Read2
	if(verbose) cat(sprintf("%4i / %4i\tREAD1>READ2\n", i, length(chimeric.combine.gr.tmp)))
	sequences     = get.fusion.sequence(donor = chimeric.combine.gr.tmp[i], acceptor = chimeric.combine.gr.tmp[i]$acc)
	sequences.all = c(sequences.all, lapply(names(sequences), function(n){ c(Id = n, Donor = as.character(sequences[[n]]$Donor), Acceptor = as.character(sequences[[n]]$Acceptor)) }))
	ids           = get.peptides(sequences = sequences, ids = ids, order = "R1R2")
}


sequences.all.rbind           = as.data.frame(do.call(rbind, sequences.all))
colnames(sequences.all.rbind) = c("Id", "Donor", "Acceptor")
ids.collapse                  = do.call(rbind, lapply(ids, paste, collapse = "/"))
ids.collapse                  = as.data.frame(ids.collapse[nchar(rownames(ids.collapse)) >= size, ])
colnames(ids.collapse)        = "Name"
ids.collapse$ID               = 1:nrow(ids.collapse)
ids.collapse$Sequence         = rownames(ids.collapse)
fusions                       = AAStringSet(ids.collapse$Sequence)
names(fusions)                = ids.collapse$ID
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

##### Saving data #####
#--------------------------------------------------------------------------------


sequences.file = paste0(prefix, "_Fusions_chim", 2, ".junc", minjunc, ".size", size, ".genomic.txt")
if(verbose) cat("Write nucleotide sequences in", sequences.file)
write.table(x = sequences.all.rbind, file = sequences.file, quote = F, sep = "\t", col.names = T, row.names = F)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")



fusions.file = paste0(prefix, "_Fusions_chim", 2, ".junc", minjunc, ".size", size, ".fasta")
if(verbose) cat("Write fusion sequences in", fusions.file)
writeXStringSet(x = fusions, filepath = fusions.file, format="fasta")
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")



ids.file = paste0(prefix, "_Fusions_chim", 2, ".junc", minjunc, ".size", size, ".ids.txt")
if(verbose) cat("Write fusion ids in", ids.file)
write.table(x = ids.collapse[, c("ID", "Name")], file = ids.file, sep = "\t", col.names = F, row.names = F, quote = F)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")


#OPTIONAL: saving all in an Rdata

rdata.file = paste0(prefix, "_Fusions_chim", 2, ".junc", minjunc, ".size", size, ".RData")
if(verbose) cat("Save data in", rdata.file)
t = Sys.time()
save(list = ls(all = T), file = rdata.file)
if(verbose) cat("\tdone in ", round(difftime(Sys.time(), t, units = 'sec'), 2), "s\n", sep = "")

if(verbose) cat("\nDone in ", round(difftime(Sys.time(), tall, units = 'sec'), 2), "s\n", sep = "")
