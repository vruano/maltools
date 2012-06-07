prefixes = c("../GATK/3D7.HB3","../BWA2SOM/HMM/3D7.HB3","Data.Type/3D7.HB3","../CORTEX/3D7.HB3");
cat.names =c("GATK","BWA2SOM","PGV","CORTEX");
#cat.names = c("Original", "OK", "Coverage", "WCSA", "Biallelic", "Uniqueness", "Typable")
#prefixes = c("Data.Orig/3D7.HB3","Data.OK/3D7.HB3","Data.Cov/3D7.HB3","Data.WCSA/3D7.HB3","Data.Bi/3D7.HB3","Data.Uniq/3D7.HB3","Data.Type/3D7.HB3");
#uniq = read.table("uniq.data",header = FALSE,sep = "\t");
#cov   = read.table("original.coverage",header =TRUE, sep ="\t");
#total.cov = rowSums(cov[,4:7]);
#total.cov = cbind(cov[,1:2],total.cov);
#rm("cov");

#cov.quant = quantile(10^total.cov[,3],seq(0.05,1,by =.05));
num.cat= length(prefixes);
colors = palette(rainbow(num.cat, s = 0.5, v = 0.75))

cov.seq = seq(0.02,1,by=0.02)
cov.quant = quantile(total.cov[,3],cov.seq,names = FALSE);

num.snp.n  = rep(0,num.cat);
num.seg.n  = rep(0,num.cat);
num.non.n  = rep(0,num.cat);
num.err.n  = rep(0,num.cat);

num.snp.c  = rep(0,num.cat);
num.seg.c  = rep(0,num.cat);
num.non.c  = rep(0,num.cat);
num.err.c  = rep(0,num.cat);

#postscript("../Figures/men.sites.uniq.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/men.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.men.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.men.sites.uniq.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);

#cnc = rep(list(),14)
for (r in 1:14)
{
#	cnc[[r]] =read.table(paste("CNC/coding.MAL",r,sep=""),header = FALSE,sep=",");
}

valid.uni = seq(18,44,by=2);
num.uniq.class = length(valid.uni);
num.cov.class = length(cov.quant);

devs = rep(0,6);

for (r in 1:num.cat)
{
	prefix = prefixes[r];

	print(r);
	## number of mendelian errors at each coverage level
	num.err.cov = rep(0,length(cov.seq));

	## number of segregating sites, non-segregating sites
	len.seg.cov = rep(0,length(cov.seq));
	len.err.cov = rep(0,length(cov.seq));

	num.err.uniq.c = rep(0,num.uniq.class);
	num.tot.uniq.c = rep(0,num.uniq.class);
	len.seg.uniq.c = rep(0,num.uniq.class);
	len.err.uniq.c = rep(0,num.uniq.class);

	num.err.uniq.n = rep(0,num.uniq.class);
	num.tot.uniq.n = rep(0,num.uniq.class);
	len.seg.uniq.n = rep(0,num.uniq.class);
	len.err.uniq.n = rep(0,num.uniq.class);

	num.err.cov.c = rep(0,num.cov.class);
	num.tot.cov.c = rep(0,num.cov.class);
	len.seg.cov.c = rep(0,num.cov.class);
	len.err.cov.c = rep(0,num.cov.class);

	num.err.cov.n = rep(0,num.cov.class);
	num.tot.cov.n = rep(0,num.cov.class);
	len.seg.cov.n = rep(0,num.cov.class);
	len.err.cov.n = rep(0,num.cov.class);

	for (j in 1:14)
	{
		print(j)
		## FILE NAMES
		if ((prefix == "../GATK/3D7.HB3")|(prefix == "../BWA2SOM/HMM/3D7.HB3"))
		{
			errant.file = paste(prefix,"errant.MAL.",sep="/");
			errant.file = paste(errant.file,j,sep="");
			errant.file = paste(errant.file,".txt",sep="");
		}
		else
		{
			errant.file = paste(prefix,"errant.MAL",sep="/");
			errant.file = paste(errant.file,j,sep="");
			errant.file = paste(errant.file,".txt",sep="");

		}
		print(errant.file);
		positions.file = paste(prefix,"positions.MAL",sep="/");
		positions.file = paste(positions.file,j,sep="");

		parents.file = paste(prefix,"parents.MAL",sep="/");
		parents.file = paste(parents.file,j,sep="");

		progeny.file = paste(prefix,"progeny.MAL",sep="/");
		progeny.file = paste(progeny.file,j,sep="");

		errants   = read.table(errant.file,header =FALSE, sep ="\t");

		positions = read.table(positions.file);

		parents   = read.table(parents.file,header = FALSE, sep =",");
		progeny   = read.table(progeny.file,header = FALSE, sep =",");

	#	print("read.");
		status = cnc[[1]][positions[,1],1];
		## TOTAL DATA, NUM SAMPLES, ROW SUMS
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
	
#		non.sites = intersect(non.sites, which( row.sum!=0 & row.sum!=num.sam ) ) 

		## CALCULATE 
		num.snp.c[r]  = num.snp.c[r]  + length(coding.sites);
		num.seg.c[r]  = num.seg.c[r]  + length(seg.sites.c);
		num.non.c[r]  = num.non.c[r]  + length(non.sites.c);
		num.err.c[r]  = num.err.c[r]  + sum(errants[coding.sites,][seg.sites.c,])

		if (length(noncod.sites)>0)
		{
			num.snp.n[r]  = num.snp.n[r]  + length(noncod.sites);
			num.seg.n[r]  = num.seg.n[r]  + length(seg.sites.n);
			num.non.n[r]  = num.non.n[r]  + length(non.sites.n);
			num.err.n[r]  = num.err.n[r]  + sum(errants[noncod.sites,][seg.sites.n,])
		}

		#print("init calcs...")

#		var.sites = intersect(which(parents[,1]!=parents[,2] ), which( row.sum!=0 & row.sum!=num.sam ) ) 
		## MENDELS, ERRANTS, POSITIONS by VAR SITES
		errants.c   = errants[coding.sites,][seg.sites.c,]
		pos.seg.c = as.matrix(positions[coding.sites,1][seg.sites.c])
		pos.err.c = as.matrix(positions[coding.sites,1][non.sites.c])
		
		if (length(noncod.sites)>0)
		{
			errants.n   = errants[noncod.sites,][seg.sites.n,]		
			pos.seg.n = as.matrix(positions[noncod.sites,1][seg.sites.n])
			pos.err.n = as.matrix(positions[noncod.sites,1][non.sites.n])
		}



		#print("more calcs...")	
		## PULL OUT CHROMO-SPECIFIC UNIQ AND COV VALUES
		chromo    = paste("MAL",j,sep ="")
		uniq.mal  = uniq[which(uniq[,1] == chromo),]
		cov.mal   = total.cov[which(total.cov[,1] == chromo),];

#		## TRANSLATE TO POSITION AND COV NUMBERING
		these.cov.seg.c = which(pos.seg.c[,1] %in% cov.mal[,2])
		those.cov.seg.c = which(cov.mal[,2] %in% pos.seg.c[,1])
		these.cov.err.c = which(pos.err.c[,1] %in% cov.mal[,2])
		those.cov.err.c = which(cov.mal[,2] %in% pos.err.c[,1])

		if (length(noncod.sites)>0)
		{
			these.cov.seg.n = which(pos.seg.n[,1] %in% cov.mal[,2])
			those.cov.seg.n = which(cov.mal[,2] %in% pos.seg.n[,1])
			these.cov.err.n = which(pos.err.n[,1] %in% cov.mal[,2])
			those.cov.err.n = which(cov.mal[,2] %in% pos.err.n[,1])
		}
		## GET THE ERR, MEN, POS and COV THAT CORRESPOND

		errant.cov.c = errants.c[these.cov.seg.c,];
		
		cov.seg.c  = cov.mal[those.cov.seg.c,];	
		errant.cov.n = errants.n[these.cov.seg.n,];

		cov.seg.n  = cov.mal[those.cov.seg.n,];
		cov.err.n  = cov.mal[those.cov.err.n,];

		#print("still more ...")	

		for (k in 1:length(cov.quant))
		{
			if (k ==1)
			{
				this.cov.err.c = which(cov.seg.c[,3] < cov.quant[k])
				this.cov.err.n = which(cov.seg.n[,3] < cov.quant[k])
			}
			else
			{
				this.cov.err.c = intersect(which(cov.seg.c[,3] < cov.quant[k]), which(cov.seg.c[,3] >= cov.quant[k-1]))
				this.cov.err.n = intersect(which(cov.seg.n[,3] < cov.quant[k]), which(cov.seg.n[,3] >= cov.quant[k-1]))
			}

			if (length(this.cov.err.c)>0)
			{
				len.err.cov.c[k] = len.err.cov.c[k] + length(this.cov.err.c);
				num.err.cov.c[k] = num.err.cov.c[k] + sum(errant.cov.c[this.cov.err.c,]);
			}
			if ((length(noncod.sites)>0)&(length(this.cov.err.n)>0))
			{
				len.err.cov.n[k] = len.err.cov.n[k] + length(this.cov.err.n);						
				num.err.cov.n[k] = num.err.cov.n[k] + sum(errant.cov.n[this.cov.err.n,]);
			}
		}

		#print("done");

		these.uniq.seg.c = which(pos.seg.c[,1] %in% uniq.mal[,2])
		those.uniq.seg.c = which(uniq.mal[,2] %in% pos.seg.c[,1])

		these.uniq.seg.n = which(pos.seg.n[,1] %in% uniq.mal[,2])
		those.uniq.seg.n = which(uniq.mal[,2] %in% pos.seg.n[,1])


		errant.uniq.c = errants.c[these.uniq.seg.c,];
		unique.err.c = uniq.mal[those.uniq.seg.c,];

		errant.uniq.n = errants.n[these.uniq.seg.n,];
		unique.err.n = uniq.mal[those.uniq.seg.n,];


	#	print("still more ...")	
		for (k in 1:length(valid.uni))
		{
			this.uniq.err.c = which(unique.err.c[,3]==valid.uni[k])
			this.uniq.err.n = which(unique.err.n[,3]==valid.uni[k])

			if (length(this.uniq.err.c)>0)
			{
				len.err.uniq.c[k] = len.err.uniq.c[k] + length(this.uniq.err.c);
				num.err.uniq.c[k] = num.err.uniq.c[k] + sum(errant.uniq.c[this.uniq.err.c,]);
			}
			if ((length(noncod.sites)>0)&(length(this.uniq.err.n)>0))
			{
				len.err.uniq.n[k] = len.err.uniq.n[k] + length(this.uniq.err.n);
				num.err.uniq.n[k] = num.err.uniq.n[k] + sum(errant.uniq.n[this.uniq.err.n,]);
			}
		}
	} 

#	print("plot...");
	if(r == 1)
	{
		postscript("../Figures/err.pipe.uniq.coding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,num.err.uniq.c/(20*len.err.uniq.c),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Errant call rate", xlab  = "Uniqueness score", ,ylim=c(0,0.17),main = "Coding")
		devs[1] = dev.cur();
		postscript("../Figures/err.pipe.uniq.noncoding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,num.err.uniq.n/(20*len.err.uniq.n),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Errant call rate", xlab  = "Uniqueness score", ,ylim=c(0,0.17), main = "Noncoding")
		devs[2] = dev.cur();
		postscript("../Figures/err.pipe.uniq.both.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,(num.err.uniq.c+num.err.uniq.n)/(20*(len.err.uniq.c+len.err.uniq.n)),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Errant call rate", xlab  = "Uniqueness score", ,ylim=c(0,0.17), main = "All")
		devs[3] = dev.cur();
		postscript("../Figures/err.pipe.cov.coding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),num.err.cov.c/(20*len.err.cov.c),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.2),xlab = "Percentile coverage",ylab="Errant call rate",axes=FALSE,main = "Coding");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[4] = dev.cur();
		postscript("../Figures/err.pipe.cov.noncoding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),num.err.cov.n/(20*len.err.cov.n),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.2),xlab = "Percentile coverage",ylab="Errant call rate",axes=FALSE, main = "Noncoding");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[5] = dev.cur();
		postscript("../Figures/err.pipe.cov.both.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),(num.err.cov.n+num.err.cov.c)/(20*(len.err.cov.n+len.err.cov.c)),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.2),xlab = "Percentile coverage",ylab="Errant call rate",axes=FALSE, main = "All");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[6] = dev.cur();
		
	}
	else
	{
	#	print("adding...");
		dev.set(devs[1]);
		points(valid.uni,num.err.uniq.c/(20*len.err.uniq.c),type = "l",col =colors[r],pch = 19,lwd=2)

		if (length(noncod.sites)>0)
		{
			dev.set(devs[2]);
			points(valid.uni,num.err.uniq.n/(20*len.err.uniq.n),type = "l",col =colors[r],pch = 19,lwd=2)
		}

		dev.set(devs[3]);
		points(valid.uni,(num.err.uniq.c+num.err.uniq.n)/(20*(len.err.uniq.c+len.err.uniq.n)),type = "l",col =colors[r],pch = 19,lwd=2)

		dev.set(devs[4]);		
		points(seq(100,1,by =-2),num.err.cov.c/(20*len.err.cov.c),type = "l",col =colors[r],pch = 19,lwd=2)
	
		if (length(noncod.sites)>0)
		{
			dev.set(devs[5]);
			points(seq(100,1,by =-2),num.err.cov.n/(20*len.err.cov.n),type = "l",col =colors[r],pch = 19,lwd=2)
		}

		dev.set(devs[6]);
		points(seq(100,1,by =-2),(num.err.cov.n+num.err.cov.c)/(20*(len.err.cov.n+len.err.cov.c)),type = "l",col =colors[r],pch = 19,lwd=2)
	}
}
#print("done.");

for (i in 1:length(devs))
{
	dev.set(devs[i]);
	legend(x ="topleft", legend = cat.names, col = colors, pch = 19,cex = 1.3)
	if ((i == 1)|(i==2)|(i==3)|(i==4))
	{

#		legend(x ="topleft", legend = c("Original","OK","Good coverage","WCSA","Biallelic","Uniqueness", "Typable"), col = colors, pch = 19,cex = 1.3)
	}
	if (i == 5)
	{
#		legend(x ="bottomright", legend = c("Original","OK","Good coverage","WCSA","Biallelic","Uniqueness", "Typable"), col = colors, pch = 19,cex = 1.3)
	}
	if (i == 6)
	{
#		legend(x ="topleft", legend = c("Original","OK","Good coverage","WCSA","Biallelic","Uniqueness", "Typable"), col = colors, pch = 19,cex = 1.2)
	}
}


#out.mat.n = as.matrix(rbind(num.snp.n,num.seg.n,num.non.n,num.err.n));
#colnames(out.mat.n) = c("SNPs", "Seg", "Non","Err");
#rownames(out.mat.n) =cat.names;
#write.table(out.mat.n,"out.mat.err.n", quote= FALSE, sep= " & ");

for (i in 1:length(devs))
{
	dev.off(devs[i])
}
