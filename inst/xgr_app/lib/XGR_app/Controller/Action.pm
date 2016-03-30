package XGR_app::Controller::Action;
use XGR_app::Controller::Utils;
use Mojo::Base 'Mojolicious::Controller';
use JSON;

# Render template "index.html.ep"
# Render template "demo.html.ep"
sub index {
	my $c = shift;
  	$c->render();
}

# Render template "enricherGenes.html.ep"
sub enricher_genes {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'DO'; # by default: DO
  	my $genelist = $c->param('genelist');
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $ontology_algorithm = $c->param('ontology_algorithm');
  	
  	my $true_path_rule='FALSE';
  	if(defined($ontology_algorithm)){
  		if($ontology_algorithm ne 'none'){
  			$true_path_rule='TRUE';	
  		}
  	}
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Genes.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $genelist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			XGR_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
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
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", ontology.algorithm=\"$ontology_algorithm\", true.path.rule=\"$true_path_rule\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			## notes: replace 'public/' with '/'
			$ajax_json_file=$output_filename;
			$ajax_json_file=~s/^public//g;
			$ajax_json_file=~s/.txt$/.json/g;
			print STDERR "JSON locates at $ajax_json_file\n";
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
	}
	
  	$c->render();

}

# Render template "enricherSNPs.html.ep"
sub enricher_snps {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'EF'; # by default: EF
	
  	my $population = $c->param('pop') || 'NA'; # by default: NA
  	my $r2 = $c->param('r2') || '0.8'; # by default: NA
  	
  	my $snplist = $c->param('snplist');
  	my $uploadfile = $c->req->upload('uploadfile');
	if(defined($uploadfile)) {
		$uploadfile = $uploadfile->slurp;
	}else{
		$uploadfile='';
	}
  	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
	my $ontology_algorithm = $c->param('ontology_algorithm');
  	
  	my $true_path_rule='TRUE';
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $population.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $snplist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $background_filename="";
		if($uploadfile ne ''){
			my $my_background;
			foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $uploadfile)) {
				next if($line=~/^\s*$/);
				$my_background.=$line."\n";
			}
			$background_filename=$tmpFolder.'/'.'background.'.$rand_file.'.txt';
			XGR_app::Controller::Utils::export_to_file($background_filename, $my_background);
		}
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", background.file="", output.file="", ontology="", population="", LD.r2="", size_range_min="", size_range_max="", min_overlap="", test="", ontology.algorithm="", true.path.rule="", RData.location="", ...){
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

	LD.r2 <- as.numeric(LD.r2)
	
	# perform enrichment analysis
	eTerm <- xEnricherSNPs(data=data, background=background, ontology=ontology, include.LD=population, LD.r2=LD.r2, size.range=size.range, min.overlap=min.overlap, test=test, ontology.algorithm=ontology.algorithm, true.path.rule=true.path.rule, RData.location=RData.location, ...)
	
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
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(input.file=\"$input_filename\", background.file=\"$background_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", population=\"$population\", LD.r2=\"$r2\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\", ontology.algorithm=\"$ontology_algorithm\", true.path.rule=\"$true_path_rule\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			## notes: replace 'public/' with '/'
			$ajax_json_file=$output_filename;
			$ajax_json_file=~s/^public//g;
			$ajax_json_file=~s/.txt$/.json/g;
			print STDERR "JSON locates at $ajax_json_file\n";
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	if(defined($snplist)){
		my @lines = split(/\r\n|\n/, $snplist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
	}
	    
  	$c->render();

}

# Render template "enricherGenes.html.ep"
sub enricher_yours {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
	my $example_entity_file='public/app/examples/Pfam.txt';
	my $example_annotation_file='public/app/examples/Pfam2GO.txt';
	
  	my $entityfile = $c->req->upload('entityfile') || '';
  	my $example_entity = $c->param('example_entity') || '';
	if($entityfile=~/Mojo::Upload/ && $example_entity ne '') {
		if(defined($example_entity)){
			$entityfile = XGR_app::Controller::Utils::read_from_file($example_entity_file);
		}
	}elsif($entityfile=~/Mojo::Upload/){
		$entityfile = $entityfile->slurp;
	}
	
  	my $annofile = $c->req->upload('annofile') || '';
  	my $example_anno = $c->param('example_anno') || '';
	if($annofile=~/Mojo::Upload/ && $example_anno ne '') {
		if(defined($example_anno)){
			$annofile=XGR_app::Controller::Utils::read_from_file($example_annotation_file);
		}
	}elsif($annofile=~/Mojo::Upload/){
		$annofile = $annofile->slurp;
  	}
  	
	my $size_range_min = $c->param('size_range_min') || 1;
	my $size_range_max = $c->param('size_range_max') || 1000000;
	my $min_overlap = $c->param('min_overlap') || 1;
	
	my $test = $c->param('test');
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
  	
  	if($entityfile ne '' and $annofile ne ''){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = 'Custom.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Yours.'.$rand_file.'.txt';
		my $anno_filename=$tmpFolder.'/'.'annotation.Yours.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'enrichment.Yours.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'enrichment.Yours.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $entityfile)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $my_anno;
		foreach my $line (split(/\n/, $annofile)) {
			$my_anno.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($anno_filename, $my_anno);
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (data.file="", annotation.file="", output.file="", size_range_min="", size_range_max="", min_overlap="", test="", ...){

	size.range <- as.numeric(c(size_range_min, size_range_max))
	min.overlap <- as.numeric(min_overlap)

	# perform enrichment analysis
	eTerm <- xEnricherYours(data.file=data.file, annotation.file=annotation.file, size.range=size.range, min.overlap=min.overlap, test=test, ...)
	
	if(class(eTerm)=="eTerm"){
		# save enrichment results to the output file
		res <- xEnrichViewer(eTerm, top_num=length(eTerm$adjp), sortBy="adjp", details=TRUE)
		res_f <- res
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
	
		# save to the json file
		res_f <- toJSON(res_f, pretty=T, digits=10)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
R_pipeline(data.file=\"$input_filename\", annotation.file=\"$anno_filename\", output.file=\"$output_filename\", size_range_min=\"$size_range_min\", size_range_max=\"$size_range_max\", min_overlap=\"$min_overlap\", test=\"$test\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			## notes: replace 'public/' with '/'
			$ajax_json_file=$output_filename;
			$ajax_json_file=~s/^public//g;
			$ajax_json_file=~s/.txt$/.json/g;
			print STDERR "JSON locates at $ajax_json_file\n";
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	    
  	$c->render();

}

# Render template "socialiserGenes.html.ep"
sub socialiser_genes {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'DO'; # by default: DO
  	my $genelist = $c->param('genelist');
  	
  	my $method_term = $c->param('method_term');
	my $measure = $c->param('measure');
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $ontology.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'similarity.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'similarity.Genes.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $genelist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", ontology="", measure="", method.term="", RData.location="", ...){
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	# perform enrichment analysis
	sim <- xSocialiserGenes(data=data, ontology=ontology, measure=measure, method.term=method.term, RData.location=RData.location, ...)
	
	if(class(sim)=="igraph"){
		# save similarity results to the output file
		res_f <- igraph::get.data.frame(sim, what="edges")
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
		
		# save to the json file (for dataTable)
		res_f <- toJSON(res_f, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# save to the json file (for highCharts)
		# reordering via hierarchical clustering
		s <- igraph::get.adjacency(sim, type="both", attr="weight", edges=F, names=T,sparse=F)
		#s <- xConverter(sim, from="igraph", to="dgCMatrix")
		data <- matrix(as.numeric(s),nrow=nrow(s))
		data[is.na(data)] <- 0
		colnames(data) <- rownames(data) <- colnames(s)
        distance <- as.dist(supraHex::sDistance(data, metric=c("pearson","spearman","kendall","euclidean","manhattan","cos","mi")[4]))
        cluster <- hclust(distance, method="average")
        ordering <- cluster$order
		new_data <- data[ordering, ordering]
		#supraHex::visHeatmapAdv(new_data, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.5), cexRow=0.7, cexCol=0.7, KeyValueName="Semantic similarity")
		max_n <- min(nrow(new_data), 30)
		m <- new_data[1:max_n, 1:max_n]
		d <- data.frame(
				i=rep(0:(nrow(m)-1),ncol(m)), 
				j=rep(0:(ncol(m)-1),each=nrow(m)),
                val=as.vector(m)
				)
		d <- as.matrix(d)
		res_f <- toJSON(list(data=d, category=colnames(m)), pretty=T)
		output.file.json <- gsub(".txt$", ".json2", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".pdf"
		output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
		pdf(file=output.file.pdf, height=8, width=8, compress=T)
		xCircos(g=sim, entity="Gene", top_num=50, ideogram=F, entity.label.cex=0.8, RData.location=RData.location)
		dev.off()
		
	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
library(RCircos)
library(GenomicRanges)
R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", measure=\"$measure\", method.term=\"$method_term\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			## notes: replace 'public/' with '/'
			$ajax_json_file=$output_filename;
			$ajax_json_file=~s/^public//g;
			$ajax_json_file=~s/.txt$/.json/g;
			print STDERR "JSON locates at $ajax_json_file\n";
			
			### for pdf
			$ajax_pdf_file=$output_filename;
			$ajax_pdf_file=~s/^public//g;
			$ajax_pdf_file=~s/.txt$/.pdf/g;
			print STDERR "Circos PDF locates at $ajax_pdf_file\n";
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $ajax_pdf_file;
	$c->stash(ajax_pdf_file => $ajax_pdf_file);
	
	if(defined($genelist)){
		my @lines = split(/\r\n|\n/, $genelist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
	}
	    
  	$c->render();

}

# Render template "socialiserSNPs.html.ep"
sub socialiser_snps {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $ontology = $c->param('ontology') || 'EF'; # by default: EF
	
  	my $population = $c->param('pop') || 'NA'; # by default: NA
  	my $r2 = $c->param('r2') || '0.8'; # by default: NA
  	
  	my $snplist = $c->param('snplist');
  	my $method_term = $c->param('method_term');
	my $measure = $c->param('measure');
  	
  	my $true_path_rule='TRUE';
  	
  	# The output json file (default: '')
	my $ajax_json_file='';
	
	# The output pdf file (default: '')
	my $ajax_pdf_file='';
	
  	if(defined($snplist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $population.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'similarity.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'similarity.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n|\s+|\,|\;/, $snplist)) {
			next if($line=~/^\s*$/);
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", ontology="", population="", LD.r2="", measure=measure, method.term=method.term, RData.location="", ...){
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1]
	
	LD.r2 <- as.numeric(LD.r2)
	
	# perform enrichment analysis
	sim <- xSocialiserSNPs(data=data, ontology=ontology, include.LD=population, LD.r2=LD.r2, measure=measure, method.term=method.term, RData.location=RData.location, ...)
	
	if(class(sim)=="igraph"){
		# save similarity results to the output file
		res_f <- igraph::get.data.frame(sim, what="edges")
		utils::write.table(res_f, file=output.file, sep="\t", row.names=FALSE)
		
		# save to the json file (for dataTable)
		res_f <- toJSON(res_f, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# save to the json file (for highCharts)
		# reordering via hierarchical clustering
		s <- igraph::get.adjacency(sim, type="both", attr="weight", edges=F, names=T,sparse=F)
		#s <- xConverter(sim, from="igraph", to="dgCMatrix")
		data <- matrix(as.numeric(s),nrow=nrow(s))
		data[is.na(data)] <- 0
		colnames(data) <- rownames(data) <- colnames(s)
        distance <- as.dist(supraHex::sDistance(data, metric=c("pearson","spearman","kendall","euclidean","manhattan","cos","mi")[4]))
        cluster <- hclust(distance, method="average")
        ordering <- cluster$order
		new_data <- data[ordering, ordering]
		#supraHex::visHeatmapAdv(new_data, Rowv=F, Colv=F, dendrogram="none", colormap="white-lightpink-darkred", zlim=c(0,1.5), cexRow=0.7, cexCol=0.7, KeyValueName="Semantic similarity")
		max_n <- min(nrow(new_data), 30)
		m <- new_data[1:max_n, 1:max_n]
		d <- data.frame(
				i=rep(0:(nrow(m)-1),ncol(m)), 
				j=rep(0:(ncol(m)-1),each=nrow(m)),
                val=as.vector(m)
				)
		d <- as.matrix(d)
		res_f <- toJSON(list(data=d, category=colnames(m)), pretty=T)
		output.file.json <- gsub(".txt$", ".json2", output.file, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".pdf"
		output.file.pdf <- gsub(".txt$", ".pdf", output.file, perl=T)
		pdf(file=output.file.pdf, height=10, width=10, compress=T)
		xCircos(g=sim, entity="SNP", top_num=50, ideogram=F, entity.label.cex=0.8, RData.location=RData.location)
		dev.off()
	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(jsonlite)
library(RCircos)
library(GenomicRanges)
R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", ontology=\"$ontology\", population=\"$population\", LD.r2=\"$r2\", measure=\"$measure\", method.term=\"$method_term\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
		if(! -e $output_filename){
			print STDERR "Cannot find $output_filename\n";
		}else{
			## notes: replace 'public/' with '/'
			$ajax_json_file=$output_filename;
			$ajax_json_file=~s/^public//g;
			$ajax_json_file=~s/.txt$/.json/g;
			print STDERR "JSON locates at $ajax_json_file\n";
			
			### for pdf
			$ajax_pdf_file=$output_filename;
			$ajax_pdf_file=~s/^public//g;
			$ajax_pdf_file=~s/.txt$/.pdf/g;
			print STDERR "Circos PDF locates at $ajax_pdf_file\n";
		}
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file;
	$c->stash(ajax_json_file => $ajax_json_file);
	
	# stash $ajax_pdf_file;
	$c->stash(ajax_pdf_file => $ajax_pdf_file);
	
	if(defined($snplist)){
		my @lines = split(/\r\n|\n/, $snplist);
		my %rec_genes;
		foreach my $line (@lines) {
			next if($line=~/^\s*$/);
			my $rec;
			$rec->{ID}=$line;
			$rec->{Name}=$line;
			$rec_genes{$line}=$rec;
		}  	
		$c->stash(rec_genes => \%rec_genes);
	}
	    
  	$c->render();

}

# Render template "subneterGenes.html.ep"
sub subneter_genes {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING
  	my $genelist = $c->param('genelist');
  	
  	#my $network_confidence = $c->param('network_confidence');
	my $subnet_significance = $c->param('subnet_significance');
	my $subnet_size = $c->param('subnet_size');
	my $flag_seeds = $c->param('flag_seeds');
  	
  	# The output json file (default: '')
	my $ajax_json_file_nodes='';
	my $ajax_json_file_edges='';
  	# The output json file (default: '')
	my $ajax_txt_file_nodes='';
	my $ajax_txt_file_edges='';
	
	# The output pdf file (default: '')
	my $ajax_pdf_file_net='';
	my $ajax_pdf_file_circos='';
  	
  	if(defined($genelist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.Genes.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'subnet.Genes.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'subnet.Genes.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $genelist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", network="", seed.genes="", subnet.significance="", subnet.size="", RData.location="", ...){
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1:2]
	
	subnet.significance <- as.numeric(subnet.significance)
	if(subnet.size=="NULL"){
		subnet.size <- NULL
	}else{
		subnet.size <- as.numeric(subnet.size)
	}
	
	if(seed.genes=="Yes"){
		seed.genes <- T
	}else{
		seed.genes <- F
	}
	
	# perform network analysis
	subnet <- xSubneterGenes(data=data, network=network, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, RData.location=RData.location, ...)
	
	if(class(subnet)=="igraph"){
		# output file ".edges.txt"
		res_edges <- igraph::get.data.frame(subnet,what="edges")
		output.file.edges <- gsub(".txt$", ".edges.txt", output.file, perl=T)
		utils::write.table(res_edges, file=output.file.edges, sep="\t", row.names=F, quote=F)
		# save to the ".edges.json.txt" file (for dataTable)
		res_f <- toJSON(res_edges, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file.edges, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".nodes.txt"
		res_nodes <- igraph::get.data.frame(subnet,what="vertices")
		output.file.nodes <- gsub(".txt$", ".nodes.txt", output.file, perl=T)
		utils::write.table(res_nodes, file=output.file.nodes, sep="\t", row.names=F, quote=F)
		# save to the ".nodes.json.txt" file (for dataTable)
		res_f <- toJSON(res_nodes, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file.nodes, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".net.pdf"
		output.file.net.pdf <- gsub(".txt$", ".net.pdf", output.file, perl=T)
		width <- height <- max(vcount(subnet)/50*7, 7)
		pdf(output.file.net.pdf, width=width, height=height, compress=T)
		## do visualisation with nodes colored according to the significance (you provide)
		pattern <- -log10(as.numeric(V(subnet)$significance))
		pattern[is.infinite(pattern)] <- max(pattern[!is.infinite(pattern)])
		vmax <- ceiling(stats::quantile(pattern, 0.8))
		vmin <- floor(min(pattern))
		glayout <- layout_(subnet, with_kk())
		xVisNet(g=subnet, pattern=pattern, glayout=glayout, vertex.shape="sphere", colormap="yr", zlim=c(vmin,vmax), edge.arrow.size=0.3, newpage=F, vertex.label.color="blue", vertex.label.dist=0.35, vertex.label.font=2)
		dev.off()
		
		# output file ".circos.pdf"
		output.file.circos.pdf <- gsub(".txt$", ".circos.pdf", output.file, perl=T)
		pdf(output.file.circos.pdf, width=10, height=10, compress=T)
		xCircos(g=subnet, entity="Gene", top_num=50, ideogram=F, entity.label.cex=0.8, RData.location=RData.location)
		dev.off()
	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(igraph)
library(jsonlite)
library(RCircos)
library(GenomicRanges)
R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", network=\"$network\", seed.genes=\"$flag_seeds\", subnet.significance=\"$subnet_significance\", subnet.size=\"$subnet_size\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
			## notes: replace 'public/' with '/'
			
			#############
			### for nodes
			$ajax_json_file_nodes=$output_filename;
			$ajax_json_file_nodes=~s/^public//g;
			$ajax_json_file_nodes=~s/.txt$/.nodes.json/g;
			print STDERR "Nodes JSON locates at $ajax_json_file_nodes\n";
			### for edges
			$ajax_json_file_edges=$output_filename;
			$ajax_json_file_edges=~s/^public//g;
			$ajax_json_file_edges=~s/.txt$/.edges.json/g;
			print STDERR "Edges JSON locates at $ajax_json_file_edges\n";
			
			#############
			### for nodes
			$ajax_txt_file_nodes=$output_filename;
			$ajax_txt_file_nodes=~s/^public//g;
			$ajax_txt_file_nodes=~s/.txt$/.nodes.txt/g;
			print STDERR "Nodes TXT locates at $ajax_txt_file_nodes\n";
			### for edges
			$ajax_txt_file_edges=$output_filename;
			$ajax_txt_file_edges=~s/^public//g;
			$ajax_txt_file_edges=~s/.txt$/.edges.txt/g;
			print STDERR "Edges TXT locates at $ajax_txt_file_edges\n";
			
			#############
			### for network pdf
			$ajax_pdf_file_net=$output_filename;
			$ajax_pdf_file_net=~s/^public//g;
			$ajax_pdf_file_net=~s/.txt$/.net.pdf/g;
			print STDERR "Network PDF locates at $ajax_pdf_file_net\n";
			### for circos pdf
			$ajax_pdf_file_circos=$output_filename;
			$ajax_pdf_file_circos=~s/^public//g;
			$ajax_pdf_file_circos=~s/.txt$/.circos.pdf/g;
			print STDERR "Network PDF locates at $ajax_pdf_file_circos\n";
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file_nodes;
	$c->stash(ajax_json_file_nodes => $ajax_json_file_nodes);
	# stash $ajax_json_file_edges;
	$c->stash(ajax_json_file_edges => $ajax_json_file_edges);
	
	# stash $ajax_txt_file_nodes;
	$c->stash(ajax_txt_file_nodes => $ajax_txt_file_nodes);
	# stash $ajax_txt_file_edges;
	$c->stash(ajax_txt_file_edges => $ajax_txt_file_edges);
	
	# stash $ajax_pdf_file;
	$c->stash(ajax_pdf_file_net => $ajax_pdf_file_net);
	$c->stash(ajax_pdf_file_circos => $ajax_pdf_file_circos);
	 
  	$c->render();

}

# Render template "subneterSNPs.html.ep"
sub subneter_snps {
  	my $c = shift;
	
	if($c->req->is_limit_exceeded){
		return $c->render(status => 400, json => { message => 'File is too big.' });
	}
	
  	my $snplist = $c->param('snplist');
  	my $population = $c->param('pop') || 'NA'; # by default: NA
	my $distance_max = $c->param('distance');
	
  	my $network = $c->param('network') || 'STRING_high'; # by default: STRING_highest
  	
	my $subnet_significance = $c->param('subnet_significance');
	my $subnet_size = $c->param('subnet_size');
	my $flag_seeds = $c->param('flag_seeds');
  	
  	# The output json file (default: '')
	my $ajax_json_file_nodes='';
	my $ajax_json_file_edges='';
  	# The output json file (default: '')
	my $ajax_txt_file_nodes='';
	my $ajax_txt_file_edges='';
	
	# The output pdf file (default: '')
	my $ajax_pdf_file_net='';
	my $ajax_pdf_file_circos='';
  	
  	if(defined($snplist)){
		my $tmpFolder = $XGR_app::Controller::Utils::tmpFolder;
		my $rand_flag = int rand 99999999;
		my $rand_file = $network.'.'.$rand_flag;
		my $input_filename=$tmpFolder.'/'.'data.SNPs.'.$rand_file.'.txt';
		my $output_filename=$tmpFolder.'/'.'subnet.SNPs.'.$rand_file.'.txt';
		my $rscript_filename=$tmpFolder.'/'.'subnet.SNPs.'.$rand_file.'.r';
	
		my $my_input;
		foreach my $line (split(/\r\n|\n/, $snplist)) {
			next if($line=~/^\s*$/);
			$line=~s/\s+/\t/;
			$my_input.=$line."\n";
		}
		XGR_app::Controller::Utils::export_to_file($input_filename, $my_input);
		
		my $RData_location;
		if(-e '/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0'){
			$RData_location="/Users/hfang/Sites/SVN/github/RDataCentre/XGR/1.0.0";
		}else{
			$RData_location="~/RDataCentre/XGR/1.0.0";
		}
		
##########################################
# BEGIN: R
##########################################
my $my_rscript='
#!/usr/bin/Rscript --vanilla
# Call R script, either using one of two following options:
# 1) R --vanilla < $rscript_file; 2) Rscript $rscript_file
';

# for generating R function
$my_rscript.='
R_pipeline <- function (input.file="", output.file="", population="", network="", distance.max="", seed.genes="", subnet.significance="", subnet.size="", RData.location="", ...){
	# read input file
	data <- read.delim(file=input.file, header=F, stringsAsFactors=F)[,1:2]
	
	distance.max <- as.numeric(distance.max)
	
	if(subnet.significance=="NULL"){
		subnet.significance <- NULL
	}else{
		subnet.significance <- as.numeric(subnet.significance)
	}
	
	if(subnet.size=="NULL"){
		subnet.size <- NULL
	}else{
		subnet.size <- as.numeric(subnet.size)
	}
	
	if(seed.genes=="Yes"){
		seed.genes <- T
	}else{
		seed.genes <- F
	}
	
	# perform network analysis
	subnet <- xSubneterSNPs(data=data, include.LD=population, network=network, seed.genes=seed.genes, subnet.significance=subnet.significance, subnet.size=subnet.size, RData.location=RData.location, ...)
	
	if(class(subnet)=="igraph"){
		# output file ".edges.txt"
		res_edges <- igraph::get.data.frame(subnet,what="edges")
		output.file.edges <- gsub(".txt$", ".edges.txt", output.file, perl=T)
		utils::write.table(res_edges, file=output.file.edges, sep="\t", row.names=F, quote=F)
		# save to the ".edges.json.txt" file (for dataTable)
		res_f <- toJSON(res_edges, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file.edges, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".nodes.txt"
		res_nodes <- igraph::get.data.frame(subnet,what="vertices")
		output.file.nodes <- gsub(".txt$", ".nodes.txt", output.file, perl=T)
		utils::write.table(res_nodes, file=output.file.nodes, sep="\t", row.names=F, quote=F)
		# save to the ".nodes.json.txt" file (for dataTable)
		res_f <- toJSON(res_nodes, pretty=T)
		res_f <- paste("{","\"data\":",res_f,"}", sep="\n")
		output.file.json <- gsub(".txt$", ".json", output.file.nodes, perl=T)
		base::write(res_f, file=output.file.json)
		
		# output file ".net.pdf"
		output.file.net.pdf <- gsub(".txt$", ".net.pdf", output.file, perl=T)
		width <- height <- max(vcount(subnet)/50*7, 7)
		pdf(output.file.net.pdf, width=width, height=height, compress=T)
		## do visualisation with nodes colored according to the significance (you provide)
		pattern <- -log10(as.numeric(V(subnet)$significance))
		pattern[is.infinite(pattern)] <- max(pattern[!is.infinite(pattern)])
		vmax <- ceiling(stats::quantile(pattern, 0.8))
		vmin <- floor(min(pattern))
		glayout <- layout_(subnet, with_kk())
		xVisNet(g=subnet, pattern=pattern, glayout=glayout, vertex.shape="sphere", colormap="yr", zlim=c(vmin,vmax), edge.arrow.size=0.3, newpage=F, vertex.label.color="blue", vertex.label.dist=0.35, vertex.label.font=2)
		dev.off()
		
		# output file ".circos.pdf"
		output.file.circos.pdf <- gsub(".txt$", ".circos.pdf", output.file, perl=T)
		pdf(output.file.circos.pdf, width=10, height=10, compress=T)
		xCircos(g=subnet, entity="Gene", top_num=50, ideogram=F, entity.label.cex=0.8, RData.location=RData.location)
		dev.off()
	}
}
';

# for calling R function
$my_rscript.="
library(XGR)
library(igraph)
library(jsonlite)
library(RCircos)
library(GenomicRanges)
R_pipeline(input.file=\"$input_filename\", output.file=\"$output_filename\", population=\"$population\", network=\"$network\", distance.max=\"$distance_max\", seed.genes=\"$flag_seeds\", subnet.significance=\"$subnet_significance\", subnet.size=\"$subnet_size\", RData.location=\"$RData_location\")
";

# for calling R function
XGR_app::Controller::Utils::export_to_file($rscript_filename, $my_rscript);
if(-e $rscript_filename){
    chmod(0755, "$rscript_filename");
    
    my $command;
    if(-e '/home/hfang/R-3.2.4/bin/Rscript'){
    	$command="/home/hfang/R-3.2.4/bin/Rscript $rscript_filename";
    }else{
    	$command="/usr/local/bin/Rscript $rscript_filename";
    }
    
    if(system($command)==1){
        print STDERR "Cannot execute: $command\n";
    }else{
			## notes: replace 'public/' with '/'
			
			#############
			### for nodes
			$ajax_json_file_nodes=$output_filename;
			$ajax_json_file_nodes=~s/^public//g;
			$ajax_json_file_nodes=~s/.txt$/.nodes.json/g;
			print STDERR "Nodes JSON locates at $ajax_json_file_nodes\n";
			### for edges
			$ajax_json_file_edges=$output_filename;
			$ajax_json_file_edges=~s/^public//g;
			$ajax_json_file_edges=~s/.txt$/.edges.json/g;
			print STDERR "Edges JSON locates at $ajax_json_file_edges\n";
			
			#############
			### for nodes
			$ajax_txt_file_nodes=$output_filename;
			$ajax_txt_file_nodes=~s/^public//g;
			$ajax_txt_file_nodes=~s/.txt$/.nodes.txt/g;
			print STDERR "Nodes TXT locates at $ajax_txt_file_nodes\n";
			### for edges
			$ajax_txt_file_edges=$output_filename;
			$ajax_txt_file_edges=~s/^public//g;
			$ajax_txt_file_edges=~s/.txt$/.edges.txt/g;
			print STDERR "Edges TXT locates at $ajax_txt_file_edges\n";
			
			#############
			### for network pdf
			$ajax_pdf_file_net=$output_filename;
			$ajax_pdf_file_net=~s/^public//g;
			$ajax_pdf_file_net=~s/.txt$/.net.pdf/g;
			print STDERR "Network PDF locates at $ajax_pdf_file_net\n";
			### for circos pdf
			$ajax_pdf_file_circos=$output_filename;
			$ajax_pdf_file_circos=~s/^public//g;
			$ajax_pdf_file_circos=~s/.txt$/.circos.pdf/g;
			print STDERR "Network PDF locates at $ajax_pdf_file_circos\n";
    }
}else{
    print STDERR "Cannot find $rscript_filename\n";
}
##########################################
# END: R
##########################################
	
	}
	
	# stash $ajax_json_file_nodes;
	$c->stash(ajax_json_file_nodes => $ajax_json_file_nodes);
	# stash $ajax_json_file_edges;
	$c->stash(ajax_json_file_edges => $ajax_json_file_edges);
	
	# stash $ajax_txt_file_nodes;
	$c->stash(ajax_txt_file_nodes => $ajax_txt_file_nodes);
	# stash $ajax_txt_file_edges;
	$c->stash(ajax_txt_file_edges => $ajax_txt_file_edges);
	
	# stash $ajax_pdf_file;
	$c->stash(ajax_pdf_file_net => $ajax_pdf_file_net);
	$c->stash(ajax_pdf_file_circos => $ajax_pdf_file_circos);
	 
  	$c->render();

}


1;
