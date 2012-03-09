package MalariaGEN::AGV::Tools::CrossReportV2;

use base 'MalariaGEN::AGV::DataTemplateTool';
use strict;
use warnings;
use vars qw(%ENV);
use File::Basename qw(dirname);
use Text::Template;
use IO::File;
use Cwd qw(realpath);
use JSON::XS;
use File::Spec::Functions qw(catfile file_name_is_absolute);
use POSIX;

our $INPUTS = {
   in_dir => { type => 'file', mandatory => 1 },
   coding_tables => { type => 'file', mandatory => 1},
   uniqueness_tables => { type => 'file', mandatory => 1},
   coverage_tables => { type => 'file', mandatory => 1},
   parents => { type => 'string', multiple => 1, mandatory => 1 },
};

our $OUTPUTS = { 
   "out_dir" => { type => 'file', mandatory => 1 }
};

sub calculate_cpu_ratio {
  my ($self,%args) = @_;
  return 1;
}

sub calculate_memory {
  return 8000;
}

sub calculate_cpu_time {
  return 35 * 60;
}

sub scripts_dir {
   my ($self) = @_;
   return realpath(catfile($ENV{PGV_HOME},'resources','cross_report','Scripts'));
}

sub interpreter {
  my ($self) = @_;
  return 'Rscript';
}

1;

__DATA__
{ use File::Spec::Functions;
  @parents = @{$J->input('parents')};
  $parent_0 = $parents[0];
  $parent_1 = $parents[1];
  $coding_tables = $J->input('coding_tables');
  $uniqueness_tables = $J->input('uniqueness_tables');
  $coverage_tables = $J->input('coverage_tables');
  $output_path = $J->output('out_dir');
  $static_path = catfile($J->wd(),"Statics");
  mkdir $static_path || die "could not create statics directory '$static_path'";
  $coverage_path = catfile($static_path,"Cov");
  $uniqueness_path = catfile($static_path,"Uniq");
  $coding_path = catfile($static_path,"Coding");

  symlink($coverage_tables,$coverage_path);
  symlink($uniqueness_tables, $uniqueness_path);
  symlink($coding_tables, $coding_path);

  $input_dir = $J->input('in_dir');
  $input_path = catfile($J->wd(),"Data");
  mkdir $input_path;
  opendir($ind_fh, $input_dir) || die "could not open input directory '$input_dir'";
  while($file = readdir $ind_fh) {
      next unless -f catfile($input_dir,$file);
      if (-l catfile($input_dir,$file)) {
        $orig = catfile($input_dir,$file);
        $dest = catfile($input_path,$file);
        `cp $orig $dest`;
        $? and die "could not copy $orig to $dest";
      }
      else {
        symlink(catfile($input_dir,$file),catfile($input_path,$file)) or die "could not create symlink '$input_path/$file' to '$input_dir/$file'";
      }
  }
 
 '' }

calc.basic.stats = function(mendels,errants,parents)	
\{
	## finds the number of segregating and non-segregating sites

	seg.sites = which(parents[,1]!=parents[,2]);
	non.sites = which(parents[,1]==parents[,2]);

	## mendelian errors and errant calls 

	men.err = if (length(non.sites) == 0) 0 else sum(mendels[non.sites,]);
	err.err = if (length(seg.sites) == 0) 0 else sum(errants[seg.sites,]);

	## num samples?
	num.samples = dim(mendels)[2];

	## matrix where they go
	errors = matrix(rep(0,6),ncol=3);
	
	## write them in
	errors[1,1:2] = c(men.err,length(non.sites));
	errors[2,1:2] = c(err.err,length(seg.sites));

	## return them
	return(errors);	
\}

calc.coding.stats = function(mendels, parents, positions, cncs)
\{
	## finds the number of segregating and non-segregating sites
	non.sites = which(parents[,1]==parents[,2]);

	## num samples?
	num.samples = dim(mendels)[2];

	## matrix where they go
	coding = matrix(rep(0,6),ncol=3);

	## calculate the coding / noncoding numbers
	pos.non.seg = positions[non.sites,1];
	
	coding.status = which(cncs[,4] == 1);
	non.coding.status = which(cncs[,4] == 0);

	coding.positions  = cncs[coding.status,2];
	non.coding.positions = cncs[non.coding.status,2];


        pos.non.seg.cod = pos.non.seg[which(pos.non.seg %in% coding.positions)];
	pos.non.seg.non.cod = pos.non.seg[which(pos.non.seg %in% non.coding.positions)];

        mendels.pos = which(positions[,1] %in% pos.non.seg.cod);
       
	men.err = if (length(mendels.pos) == 0) 0 else sum(mendels[mendels.pos,]);
	coding[1,1:2] = c(men.err,length(non.sites[pos.non.seg.cod]));

        mendels.pos = which(positions[,1] %in% pos.non.seg.non.cod);
        
	men.err = if (length(mendels.pos) == 0) 0 else sum(mendels[mendels.pos,]);
	coding[2,1:2] = c(men.err,length(non.sites[pos.non.seg.non.cod]));


	## return them
	return(coding);	
\}

calc.uniq.stats = function(mendels,parents, positions, uniq)
\{

	## finds the non-segregating sites
	non.sites = which(parents[,1]==parents[,2]);
	pos.non.seg = positions[non.sites,1];

	## num samples?
	num.samples = dim(mendels)[2];

	## matrix where they go
	uniq.errors = matrix(rep(0,14*3),ncol = 3);

	## do the calculations
	uniq.seq = seq(18,44,by = 2);
	
	for (i in 1:length(uniq.seq))
	\{
		this.score = which(uniq[,3] == uniq.seq[i])
		this.score.positions = uniq[this.score,2];

		uniq.pos = intersect(pos.non.seg,this.score.positions);
                which.pos = which(positions[,1] %in% uniq.pos);
                if (length(which.pos) == 0) \{
                   uniq.errors[i,1] = 0;
                \} 
                else \{
                   uniq.errors[i,1] = sum(mendels[which.pos,]);
                \}
                uniq.errors[i,2] = length(which.pos);
	\}

	## let them go
	return(uniq.errors);
\}

calc.cov.stats = function(mendels,parents, positions, cov, quants)
\{

	## finds the non-segregating sites
	non.sites = which(parents[,1]==parents[,2]);
	pos.non.seg = positions[non.sites,1];

	## num samples?
	num.samples = dim(mendels)[2];

	## quantile vector
	cov.seq = as.vector(quants[,2]);

	## matrix where they go
	cov.errors = matrix(rep(0,length(cov.seq)*3),ncol = 3);

	## do the calculations
	for (i in 1:length(cov.seq))
	\{
		if ( i > 1)
		\{	this.score = which((cov[,3] <= cov.seq[i])&(cov[,3]>cov.seq[i-1]))	\}
		else
		\{	this.score = which(cov[,3] <= cov.seq[i]);		\}

		this.score.positions = cov[this.score,2];

		cov.pos = intersect(pos.non.seg,this.score.positions);
                which.pos = which(positions[,1] %in% cov.pos);
                if (length(which.pos) == 0) \{
                   cov.errors[i,1] = 0;
                \}
                else \{ 
		   cov.errors[i,1] = sum(mendels[which.pos,]);
		\}
                cov.errors[i,2] = length(which.pos);
	\}

	## let them go

	return(cov.errors);
\}

calc.binom.stats = function(progeny,parents)
\{
	## calculate segregating sites
	seg.sites = which(parents[,1] != parents[,2]);

	## number of samples
	num.samples = dim(progeny)[2];

	## output vector
	counts = seq(0,num.samples,by = 1);

	## calculate binomial sites
	table.out = table(rowSums(progeny[seg.sites,]==parents[seg.sites,2]));
	counts[as.numeric(names(table.out))+1] = as.vector(table.out);

	return(counts);
\}

calc.site.stats = function(parents,progeny)
\{
	## calculate number of unfixed sites

	num.unfixed = dim(progeny)[1];

	total = cbind(progeny,parents);

	row.sums = rowSums(total);

	num.ind = dim(total)[2];

	num.constant = length(which((row.sums==0)|(row.sums==num.ind)));
	num.variable = length(which((row.sums!=0)&(row.sums!=num.ind)));

	return(matrix(c(num.unfixed,num.constant,num.variable),ncol=3));

\}

plot.figures = function(out.path, figure.path)
\{

	## UNIQUENESS
	## path in
	uniq.out.path = paste(out.path,"/uniq.tab",sep="");
	uniqs = read.table(uniq.out.path,header = TRUE,sep = "&");

	## path out to figures
	uniq.fig.path = paste(figure.path,"/uniq.ps",sep="");
	postscript(uniq.fig.path, paper = "special", width = 7, height = 7,horizontal = FALSE);
	##plot
	plot(as.numeric(rownames(uniqs)), uniqs[,3],col="blue",type ="h",lwd = 2, xlab = "Nonuniqueness score",ylab ="Mendelian error rate", ylim = c(0,0.04));
	dev.off();

	## COVERAGE
	## path in
	cov.out.path = paste(out.path,"/cov.tab",sep="");
	cov = read.table(file = cov.out.path, header = TRUE, sep = "&");

	## path out to figures
	cov.fig.path = paste(figure.path,"/coverage.ps",sep="");
	postscript(cov.fig.path, paper = "special", width = 7, height = 7,horizontal = FALSE);
	##plot
	plot(as.numeric(rownames(cov)), cov[,3],col="blue",type ="h",lwd = 2, xlab = "Coverage percentile",ylab ="Mendelian error rate", ylim = c(0,0.045));
	dev.off();

	## BINOMIAL 

	## path in
	binom.out.path = paste(out.path,"/binom.tab",sep="");
	binom = read.table(file = binom.out.path, header = FALSE, sep = "&");
	## path out to figures
	binom.fig.path = paste(figure.path,"/binom.ps",sep="");
	postscript(binom.fig.path, paper = "special", width = 7, height = 7,horizontal = FALSE);

	## plot 
	z = as.numeric(unlist(binom[1,])); 
	y = as.numeric(unlist(binom[2,]));
	y.s = y/sum(y);
	max.z = max(z);
	dense = dbinom(z,max.z,0.5)
	plot(seq(0,max.z,1),dense,type="l",col="black",xlab="Counts", ylab="Density");
	polygon(x = c(seq(0,max.z,by=1),seq(max.z,0,by=-1)),y = c(dense,rep(0,length(dense))), col=grey(0.5));
	points(z, y.s,col="red",type ="h",lwd = 2, xlab = "Counts",ylab ="Mendelian error rate", ylim = c(0,0.025));
	dev.off();

\}

write.summary.stats = function(data.path, static.path, out.path, figures.path)
\{
	## output matrices
	errors = matrix(rep(0,6),ncol = 3);
	colnames(errors) =c("Errors","Sites","Rate");	
	rownames(errors) =c("Mendelian","Errant");

	coding = matrix(rep(0,6),ncol = 3);
	colnames(coding) =c("Errors","Sites","Rate");	
	rownames(coding) =c("Coding","Noncoding");

	uniqs  = matrix(rep(0,14*3),ncol = 3);
	colnames(uniqs) =c("Errors","Sites","Rate");	
	rownames(uniqs) = as.character(seq(18, 44, by=2));

	covs  = matrix(rep(0,50*3),ncol = 3);
	colnames(covs) =c("Errors","Sites","Rate");	
	rownames(covs) = as.character(seq(1,100, by=2));
	
        counts = matrix(rep(0,3),ncol = 3)
	colnames(counts) = c("Num.unfixed", "Num.constant", "Num.variable");


	## coverage quantile is same for all chromosomes
	quant.path = paste(static.path,"Cov/cov.quants",sep="/");		
	quants     = read.table(quant.path,header = FALSE, sep="\t");

	for (i in 1:14)
	\{

		## read in the chromosome specific static data
		## assumes data.path is path.to/
		## coding file is path.to/cnc.MAL1.tab
		coding.path = paste(static.path,"Coding",sep="/");		
		coding.path = paste(coding.path,"cnc.MAL",sep="/");
		coding.path = paste(coding.path,i,sep="");
		coding.path = paste(coding.path,".tab",sep="");

		cncs = read.table(coding.path,header = FALSE,sep="\t");	

		uniq.path = paste(static.path,"Uniq",sep="/");		
		uniq.path = paste(uniq.path,"uniq.data.MAL",sep="/");
		uniq.path = paste(uniq.path,i,sep="");
	
		uniq = read.table(uniq.path,header = FALSE,sep="\t");

		cov.path = paste(static.path,"Cov",sep="/");		
		cov.path = paste(cov.path,"coverage.MAL",sep="/");
		cov.path = paste(cov.path,i,sep="");
	
		cov  = read.table(cov.path, header= TRUE, sep = "\t")

		## data file paths 
		mendel.path = paste(data.path,"mendels.MAL",sep ="/");
		mendel.path = paste(mendel.path,i,sep="");
		mendel.path = paste(mendel.path,".txt",sep="");

		errant.path = paste(data.path,"errant.MAL",sep="/");
		errant.path = paste(errant.path,i,sep="");
		errant.path = paste(errant.path,"txt",sep=".");

		positions.path = paste(data.path,"positions.MAL",sep="/");
		positions.path = paste(positions.path,i,sep="");

		parents.path   = paste(data.path,"parents.MAL",sep="/");
		parents.path   = paste(parents.path,i,sep="");

		progeny.path   = paste(data.path,"progeny.MAL",sep="/");
		progeny.path   = paste(progeny.path,i,sep="");
		
		## read in all the basic data
		mendels   = read.table(mendel.path,header= FALSE,sep="\t");
		errants   = read.table(errant.path,header= FALSE,sep="\t");	
		positions = read.table(positions.path,header= FALSE,sep="\t");	
		parents   = read.table(parents.path,header= FALSE,sep=",");	
		progeny   = read.table(progeny.path,header =FALSE,sep=",");

		## small exeception: binom out vector here because it needs the number of samples
		if (i == 1)
		\{
			num.samples = dim(progeny)[2];
			binom  = matrix(rep(0,num.samples+1),ncol = (num.samples+1));
			colnames(binom) =as.character(seq(0,num.samples, by=1));
			rownames(binom) = c("Counts");
		\}

		## do the calculations
		errors = errors + calc.basic.stats(mendels,errants,parents);
		coding = coding + calc.coding.stats(mendels,parents,positions,cncs);
		uniqs  = uniqs  + calc.uniq.stats(mendels,parents,positions,uniq);
		covs   = covs   + calc.cov.stats(mendels,parents,positions,cov,quants);
		binom  = binom  + calc.binom.stats(progeny,parents);
		counts = counts + calc.site.stats(progeny,parents);

	\}

	num.samples = dim(mendels)[2];
	
        errors[1,3] = errors[1,1]/(num.samples*errors[1,2]);
        errors[2,3] = errors[2,1]/(num.samples*errors[2,2]);

        coding[1,3] = coding[1,1]/(num.samples*coding[1,2]);
        coding[2,3] = coding[2,1]/(num.samples*coding[2,2]);

        num.uniq = dim(uniqs)[1];

        for (i in 1:num.uniq)
        \{
           uniqs[i,3] = uniqs[i,1]/(num.samples*uniqs[i,2]);
        \}

        num.covs = dim(covs)[1];

        for (i in 1:num.covs)
        \{
           covs[i,3] = covs[i,1]/(num.samples*covs[i,2]);
        \}

	## output the calculations in tables
	counts.out.path = paste(out.path,"/counts.tab",sep="");
	write.table(counts,file = counts.out.path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " & ");

	error.out.path  = paste(out.path,"/errors.tab",sep="");
	write.table(errors,file = error.out.path, row.names= TRUE, col.names = TRUE, quote = FALSE, sep=" & ");
	
	coding.out.path = paste(out.path,"/coding.tab",sep="");
	write.table(coding,file = coding.out.path, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = " & ");
	
	uniq.out.path = paste(out.path,"/uniq.tab",sep="");
	write.table(uniqs,file = uniq.out.path, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = " & ");
	
	cov.out.path = paste(out.path,"/cov.tab",sep="");
	write.table(covs,file = cov.out.path, row.names = TRUE, col.names = TRUE, quote = FALSE, sep = " & ");

	binom.out.path = paste(out.path,"/binom.tab",sep="");
	write.table(binom,file = binom.out.path, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " & ");

        plot.figures(out.path,figures.path);
\}


calc.errors = function(prefix,figure,cross_1,cross_2) \{
  
  prefix.parents = paste(prefix,"parents.MAL",sep="/");
  prefix.progeny = paste(prefix,"progeny.MAL",sep="/");
  prefix.paths = paste(prefix,"paths.MAL",sep="/");
  prefix.positions = paste(prefix,"positions.MAL",sep="/");
  
  for (i in 1:14)
    \{
      mendel = paste(prefix,paste("mendels.MAL",i,".txt",sep=""),sep="/");
      errant = paste(prefix,paste("errant.MAL",i,".txt",sep=""),sep="/");
      file.name = paste(figure,paste("cross.MAL",i,".ps",sep=""),sep = "/");
      
      parent.file = paste(prefix.parents,i,sep="");
      
      print(paste("Readding parents... ",parent.file));	
      parents = read.table(parent.file,header = FALSE,sep =",");
      
      progeny.file = paste(prefix.progeny,i,sep="");
      progeny = read.table(progeny.file, header = FALSE, sep = ",");
      
      num.samples = dim(progeny)[2];
      
      paths.file = paste(prefix.paths,i,sep="");
      print(paths.file);
      paths = read.table(paths.file,header = FALSE,sep=",");
      
      positions.file = paste(prefix.positions,i,sep="");
      print(positions.file);
      positions = read.table(positions.file,header = FALSE,sep=",");
      positions = positions[,1];
      
      postscript(file.name,paper="special", width = 14, height = 10,horizontal = FALSE);
      layout(matrix(c(1,1,2,2),byrow = TRUE, ncol = 2),heights = c(3,1));
      par(mar=c(0,7,1.5,2),oma=c(0,0,2,0));
      out.mat = cbind(t(paths),(parents));
      image(as.matrix(out.mat),col = grey((0:10)/10),axes= FALSE);
      
      mendel.errors = matrix(rep(0,num.samples*length(positions)),ncol = length(positions));
      
      joint.snps = which(parents[,2] != parents[,1]);
      
      for (k in 1:num.samples)
	\{
	  mendel.errors[k,which((progeny[,k] != parents[,1])&(progeny[,k] != parents[,2]))] = 1;
	\}
      
      
      for (k in 1:num.samples)
	\{
	  num.mendel = sum(mendel.errors[k,])
	    if (num.mendel > 0)
	      \{
		x = which(mendel.errors[k,]==1)/length(mendel.errors[k,]);
		y = rep(0.002+(k-1)/(num.samples+1),num.mendel);
		
		points( x, y,col ="red",pch = 19, cex=0.5)
	      \}
	\}
   
      for(k in 0:num.samples)
        \{
          if (k == num.samples)
            \{
               lines(x= c(0,1), y = c(k/(num.samples+1)-1/((num.samples+2)*2),k/(num.samples+1)-1/((num.samples+2)*2)),col = "purple",lwd =5,cex = 1);
            \}
          else
            \{
               lines(x= c(0,1), y = c(k/(num.samples+1)-1/((num.samples+2)*2),k/(num.samples+1)-1/((num.samples+2)*2)),col = grey(0.5),lwd =1,cex = 1);
            \}
        \}

   num.positions = length(positions);
   hmm.errors = matrix(rep(0,num.samples*length(positions)),ncol = length(positions));
   
   for (k in 1:num.samples)
     \{
       p.1 = rep(0,num.positions);
       noughts = which(paths[k,]==0);
       calls.hmm = union(which(parents[noughts,1]!= progeny[noughts,k]),which(parents[-noughts,2]!=progeny[-noughts,k]))
       hmm.errors[k,calls.hmm] = 1;
       hmm.errors[k,-joint.snps]=0;
       
       if (length(calls.hmm)>0)
	 \{
	   points( x = calls.hmm/length(calls.hmm), y = rep(0.009+(k-1)/(num.samples+1),length(calls.hmm)),col ="blue",pch = 19, cex =0.5);
	 \}	
     \}
   
   print(paste("Wrinting errants ",errant));
   write.table(t(hmm.errors), file=errant, col.names = FALSE, row.names = FALSE, sep = "\t");
   print(paste("Wrinting mendel ",mendel));
   write.table(t(mendel.errors), file=mendel, col.names = FALSE, row.names = FALSE, sep = "\t");
   
   total.num = num.samples + 2;
   axis(2, labels = c(seq(1,num.samples,by = 1),cross_1,cross_2), at = seq(0.5,total.num+0.5,length.out= (total.num))/total.num, cex.axis = 0.8);
   box();
   string = paste("HMM path, with Mendelian errors and errant SNP calls: Chromosome ",i,sep = "");
   title(main = string, lty = 3);
   
   snp.positions = positions;
   num.snps = length(snp.positions);
   
   max.val = max(snp.positions);
   min.val = min(snp.positions);
   plot(seq(0,1,length.out = max.val),rep(1,max.val), type="n", axes=F, xlab= "SNP Position", ylab = "",xlim = c(0.035,0.96), ylim = c(0.03,1));
   par(mar=c(4,7,2,2),oma=c(1,0,0,0));
   box();
   
   for (k in seq(1,num.snps,length.out=min(200,length(snp.positions))))
   \{
     lines(x = c((snp.positions[k])/max.val,(k-1)/num.snps),y=c(0,1),col ="black");
   \}

   axis(1,at = c(0.04,0.25,0.5,0.75,0.96), label = round(c(0.04,0.25,0.5,0.75,0.96)*max.val));
   title( xlab = "SNP position");
   
   dev.off();
  \}
\}

print (system("echo ${YYY}",intern=T));

if (system("hmm {$input_path}/ 0.00001 0.0001") != 0) \{
  write("Failure executing hmm, command: hmm {$input_path}/ 0.00001 0.0001",stderr());
  quit(save="no",status=1);   
\}
if (system("cp {$input_path}/paths.* {$output_path}") != 0) \{
  write("Failure coping path files to output, command: cp {$input_path}/paths.* {$output_path}",stderr());
  quit(save="no",status=2);   
\}

system("mkdir -p {$output_path}/Figures");
calc.errors("{$input_path}","{$output_path}/Figures","{$parent_0}","{$parent_1}");
if (system("cp {$input_path}/mendels.* {$input_path}/errant.* {$output_path}") != 0) \{
  write("Failure coping mendels files to output, command: cp {$input_path}/mendels.* {$output_path}",stderr());
  quit(save="no",status=2);   
\}


write.summary.stats("{$input_path}","{$static_path}","{$output_path}","{$output_path}/Figures");
