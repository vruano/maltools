cat("-- reading arguments\n", sep = "");

cmd_args = commandArgs();

in_path  = cmd_args[4];
out_path = cmd_args[5];

uniq = read.table("Invariants/uniq.data",header = FALSE,sep = "\t");

num.snp.n  = 0;
num.seg.n  = 0;
num.non.n  = 0;
num.err.n  = 0;
num.men.n  = 0;

num.snp.c  = 0;
num.seg.c  = 0;
num.non.c  = 0;
num.err.c  = 0;
num.men.c  = 0;


cnc = rep(list(),14)

for (r in 1:14)
{
	cnc[[r]] =read.table(paste("Invariants/Coding/coding.MAL",r,sep=""),header = FALSE,sep=",");
}

valid.uni = seq(18,44,by=2);
num.uniq.class = length(valid.uni);

devs = rep(0,6);

for (r in 1:num.cat)
{
	prefix = prefixes[r];

	print(r);

	## number of mendelian errors at each coverage level
	num.err.uniq.c = rep(0,num.uniq.class);
	num.men.uniq.c = rep(0,num.uniq.class);

	num.tot.uniq.c = rep(0,num.uniq.class);
	len.seg.uniq.c = rep(0,num.uniq.class);
	len.non.uniq.c = rep(0,num.uniq.class);

	len.err.uniq.c = rep(0,num.uniq.class);
	len.men.uniq.c = rep(0,num.uniq.class);

	num.err.uniq.n = rep(0,num.uniq.class);
	num.men.uniq.n = rep(0,num.uniq.class);

	num.tot.uniq.n = rep(0,num.uniq.class);
	len.seg.uniq.n = rep(0,num.uniq.class);
	len.non.uniq.n = rep(0,num.uniq.class);

	len.err.uniq.n = rep(0,num.uniq.class);
	len.men.uniq.n = rep(0,num.uniq.class);

	
	for (j in 1:14)
	{
		errant.file = paste(prefix,"errant.MAL",sep="");
		errant.file = paste(errant.file,j,sep="");
		errant.file = paste(errant.file,".txt",sep="");

		mendel.file = paste(prefix,"mendels.MAL",sep="");
		mendel.file = paste(mendel.file,j,sep="");
		mendel.file = paste(mendel.file,".txt",sep="");

		positions.file = paste(prefix,"positions.MAL",sep="/");
		positions.file = paste(positions.file,j,sep="");

		parents.file = paste(prefix,"parents.MAL",sep="/");
		parents.file = paste(parents.file,j,sep="");

		progeny.file = paste(prefix,"progeny.MAL",sep="/");
		progeny.file = paste(progeny.file,j,sep="");

		errants   = read.table(errant.file,header =FALSE, sep ="\t");
		mendels   = read.table(mendel.file,header =FALSE, sep ="\t");

		positions = read.table(positions.file);
		parents   = read.table(parents.file,header = FALSE, sep =",");
		progeny   = read.table(progeny.file,header = FALSE, sep =",");
	
	}
		
}


