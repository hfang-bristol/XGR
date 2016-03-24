
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file

R_pipeline <- function (input.file="", background.file="", output.file="", ontology="", size_range_min="", size_range_max="", min_overlap="", test="", ontology.algorithm="", true.path.rule="", RData.location="", ...){
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	if(background.file!=""){
		# read background file
		background <- read.delim(file=background.file, header=F, stringsAsFactors=F)[,1]	
	}else{
		background <- NULL
	}
	
	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)
	
	# perform enrichment analysis
	eTerm <- xEnricherGenes(data=data, background=background, ontology=ontology, size.range=size.range, min.overlap=min.overlap, test=test, ontology.algorithm=ontology.algorithm, true.path.rule=true.path.rule, RData.location=RData.location, ...)
	
	if(class(eTerm)=="eTerm"){
		# save enrichment results to the output file
		res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res_f <- data.frame(term=rownames(res), res)
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
	}
}

library(XGR)
library(jsonlite)
R_pipeline(input.file="public/tmp/data.Genes.MsigdbC2CPall.41919197.txt", background.file="", output.file="public/tmp/enrichment.Genes.MsigdbC2CPall.41919197.txt", ontology="MsigdbC2CPall", size_range_min="5", size_range_max="2000", min_overlap="3", test="hypergeo", ontology.algorithm="none", true.path.rule="FALSE", RData.location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0")
