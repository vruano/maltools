cat("-- reading arguments\n", sep = "");

cmd_args = commandArgs();

in_path  = cmd_args[4];
out_path = cmd_args[5];
uq_file = cmd_args[6]; 
coding.tables = cmd_args[7];

print(paste("Reading uniqueness ",uq_file))
uniq = read.table(uq_file,header = FALSE,sep = "\t");

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
	cnc[[r]] =read.table(paste(coding.tables,"/coding.MAL",r,sep=""),header = FALSE,sep=",");
}

valid.uni = seq(18,44,by=2);
num.uniq.class = length(valid.uni);

devs = rep(0,6);


## number of mendelian errors at each uniqueness level
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

	errant.file = paste(in_path,"errant.MAL",sep="");
	errant.file = paste(errant.file,j,sep="");
	errant.file = paste(errant.file,".txt",sep="");

	mendel.file = paste(in_path,"mendels.MAL",sep="");
	mendel.file = paste(mendel.file,j,sep="");
	mendel.file = paste(mendel.file,".txt",sep="");

	positions.file = paste(in_path,"positions.MAL",sep="/");
	positions.file = paste(positions.file,j,sep="");

	parents.file = paste(in_path,"parents.MAL",sep="/");
	parents.file = paste(parents.file,j,sep="");

	progeny.file = paste(in_path,"progeny.MAL",sep="/");
	progeny.file = paste(progeny.file,j,sep="");

	errants   = read.table(errant.file,header =FALSE, sep ="\t");
	mendels   = read.table(mendel.file,header =FALSE, sep ="\t");

	positions = read.table(positions.file);
	parents   = read.table(parents.file,header = FALSE, sep =",");
	progeny   = read.table(progeny.file,header = FALSE, sep =",");

	
	status = cnc[[1]][positions[,1],1];
	
	totals    = cbind(parents,progeny);
	num.sam   = dim(totals)[2];
	row.sum   = rowSums(totals); 
		
	coding.sites = which(status == 1);
	noncod.sites = which(status == 0);

	## FIND VARIABLE SITES 
	seg.sites.c = which(parents[coding.sites,1] != parents[coding.sites,2]);
	non.sites.c = which(parents[coding.sites,1] == parents[coding.sites,2]);

	seg.sites.n = which(parents[noncod.sites,1] != parents[noncod.sites,2]);
	non.sites.n = which(parents[noncod.sites,1] == parents[noncod.sites,2]);
	

	## CALCULATE NUMBERS

	num.snp.c  = num.snp.c  + length(coding.sites);

	num.seg.c  = num.seg.c  + length(seg.sites.c);
	num.non.c  = num.non.c  + length(non.sites.c);

	num.err.c  = num.err.c  + sum(errants[coding.sites,][seg.sites.c,])
	num.men.c  = num.err.c  + sum(mendels[coding.sites,][non.sites.c,])

	if (length(noncod.sites)>0)
	{
		num.snp.n  = num.snp.n  + length(noncod.sites);

		num.seg.n  = num.seg.n  + length(seg.sites.n);
		num.non.n  = num.non.n  + length(non.sites.n);
		num.err.n  = num.err.n  + sum(errants[noncod.sites,][seg.sites.n,])
		num.men.n  = num.men.n  + sum(mendels[noncod.sites,][non.sites.n,])
	}
	
	errants.c   = errants[coding.sites,][seg.sites.c,]
	mendels.c   = mendels[coding.sites,][non.sites.c,]
	pos.seg.c = as.matrix(positions[coding.sites,1][seg.sites.c])
	pos.non.c = as.matrix(positions[coding.sites,1][non.sites.c])
		
	if (length(noncod.sites)>0)
	{
		errants.n   = errants[noncod.sites,][seg.sites.n,]	
		mendels.n   = mendels[noncod.sites,][non.sites.n,]		
		pos.seg.n = as.matrix(positions[noncod.sites,1][seg.sites.n])
		pos.non.n = as.matrix(positions[noncod.sites,1][non.sites.n])
	}

	## PULL OUT CHROMO-SPECIFIC UNIQ VALUES
	chromo    = paste("MAL",j,sep ="")
	uniq.mal  = uniq[which(uniq[,1] == chromo),]

	these.uniq.seg.c = which(pos.seg.c[,1] %in% uniq.mal[,2])
	those.uniq.seg.c = which(uniq.mal[,2] %in% pos.seg.c[,1])

	these.uniq.non.c = which(pos.non.c[,1] %in% uniq.mal[,2])
	those.uniq.non.c = which(uniq.mal[,2] %in% pos.non.c[,1])

	unique.seg.c = uniq.mal[those.uniq.seg.c,];
	unique.non.c = uniq.mal[those.uniq.non.c,];

	errant.uniq.c = errants.c[these.uniq.seg.c,];
	mendel.uniq.c = mendels.c[these.uniq.non.c,];

	print(length(noncod.sites));

	if (length(noncod.sites)>0)
	{
		these.uniq.seg.n = which(pos.seg.n[,1] %in% uniq.mal[,2])
		those.uniq.seg.n = which(uniq.mal[,2] %in% pos.seg.n[,1])

		these.uniq.non.n = which(pos.non.n[,1] %in% uniq.mal[,2])
		those.uniq.non.n = which(uniq.mal[,2] %in% pos.non.n[,1])	
	}

	if (length(noncod.sites)>0)
	{
		errant.uniq.n = errants.n[these.uniq.seg.n,];
		mendel.uniq.n = mendels.n[these.uniq.non.n,];

		unique.seg.n = uniq.mal[those.uniq.seg.n,];
		unique.non.n = uniq.mal[those.uniq.non.n,];
	}

	#print("still more ...")	

	for (k in 1:length(valid.uni))
	{
		this.uniq.seg.c = which(unique.seg.c[,3]==valid.uni[k])
		this.uniq.non.c = which(unique.non.c[,3]==valid.uni[k])

		this.uniq.seg.n = integer(0);
		this.uniq.non.n = integer(0);

		if ( length(noncod.sites)>0)
		{
			this.uniq.seg.n = which(unique.seg.n[,3]==valid.uni[k])
			this.uniq.non.n = which(unique.non.n[,3]==valid.uni[k])		
		}

		if (length(this.uniq.seg.c)>0)
		{
			len.err.uniq.c[k] = len.err.uniq.c[k] + length(this.uniq.seg.c);
			num.err.uniq.c[k] = num.err.uniq.c[k] + sum(errant.uniq.c[this.uniq.seg.c,]);
		}

		if (length(this.uniq.non.c)>0)
		{
			len.men.uniq.c[k] = len.men.uniq.c[k] + length(this.uniq.non.c);
			num.men.uniq.c[k] = num.men.uniq.c[k] + sum(mendel.uniq.c[this.uniq.non.c,]);
		}

		if ((length(noncod.sites)>0))
		{
			len.err.uniq.n[k] = len.err.uniq.n[k] + length(this.uniq.seg.n);
			if (length(this.uniq.seg.n)>0)
			{
				num.err.uniq.n[k] = num.err.uniq.n[k] + sum(errant.uniq.n[this.uniq.seg.n,]);
			}
			
			len.men.uniq.n[k] = len.err.uniq.n[k] + length(this.uniq.non.n);
			if (length(this.uniq.non.n)>0)
			{
				num.men.uniq.n[k] = num.men.uniq.n[k] + sum(mendel.uniq.n[this.uniq.non.n,]);
			}
		}

	}
}	

print("done");	

out.mat = cbind(len.err.uniq.c, len.men.uniq.c, len.err.uniq.n, len.men.uniq.n, num.err.uniq.c, num.men.uniq.c, num.err.uniq.n, num.men.uniq.n)

colnames(out.mat) = c("Len.Err.C","Len.Men.C","Len.Err.N","Len.Men.N","Num.Err.C","Num.Men.C","Num.Err.N","Num.Men.N");
rownames(out.mat) = c(valid.uni);
out_path = paste(out_path,"/out.matrix",sep="");
write.table(out.mat,out_path, , quote=FALSE, sep ="\t");
