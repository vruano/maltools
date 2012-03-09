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
num.men.n  = rep(0,num.cat);

num.snp.c  = rep(0,num.cat);
num.seg.c  = rep(0,num.cat);
num.non.c  = rep(0,num.cat);
num.men.c  = rep(0,num.cat);

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

for (r in 1:4)
{
	prefix = prefixes[r];

	print(r);
	## number of mendelian errors at each coverage level
	num.men.cov = rep(0,length(cov.seq));

	## number of segregating sites, non-segregating sites
	len.seg.cov = rep(0,length(cov.seq));
	len.men.cov = rep(0,length(cov.seq));

	num.men.uniq.c = rep(0,num.uniq.class);
	num.tot.uniq.c = rep(0,num.uniq.class);
	len.seg.uniq.c = rep(0,num.uniq.class);
	len.men.uniq.c = rep(0,num.uniq.class);

	num.men.uniq.n = rep(0,num.uniq.class);
	num.tot.uniq.n = rep(0,num.uniq.class);
	len.seg.uniq.n = rep(0,num.uniq.class);
	len.men.uniq.n = rep(0,num.uniq.class);

	num.men.cov.c = rep(0,num.cov.class);
	num.tot.cov.c = rep(0,num.cov.class);
	len.seg.cov.c = rep(0,num.cov.class);
	len.men.cov.c = rep(0,num.cov.class);

	num.men.cov.n = rep(0,num.cov.class);
	num.tot.cov.n = rep(0,num.cov.class);
	len.seg.cov.n = rep(0,num.cov.class);
	len.men.cov.n = rep(0,num.cov.class);

	for (j in 1:14)
	{
		print(j)
		## FILE NAMES
		if ((prefix == "../GATK/3D7.HB3")|(prefix == "../BWA2SOM/HMM/3D7.HB3"))
		{
			mendel.file = paste(prefix,"mendels.MAL.",sep="/");
			mendel.file = paste(mendel.file,j,sep="");
			mendel.file = paste(mendel.file,".txt",sep="");
		}
		else
		{
			mendel.file = paste(prefix,"mendels.MAL",sep="/");
			mendel.file = paste(mendel.file,j,sep="");
			mendel.file = paste(mendel.file,".txt",sep="");

		}
		print(mendel.file);
		positions.file = paste(prefix,"positions.MAL",sep="/");
		positions.file = paste(positions.file,j,sep="");

		parents.file = paste(prefix,"parents.MAL",sep="/");
		parents.file = paste(parents.file,j,sep="");

		progeny.file = paste(prefix,"progeny.MAL",sep="/");
		progeny.file = paste(progeny.file,j,sep="");

		mendels   = read.table(mendel.file,header =FALSE, sep ="\t");

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
		num.men.c[r]  = num.men.c[r]  + sum(mendels[coding.sites,][non.sites.c,])

		if (length(noncod.sites)>0)
		{
			num.snp.n[r]  = num.snp.n[r]  + length(noncod.sites);
			num.seg.n[r]  = num.seg.n[r]  + length(seg.sites.n);
			num.non.n[r]  = num.non.n[r]  + length(non.sites.n);
			num.men.n[r]  = num.men.n[r]  + sum(mendels[noncod.sites,][non.sites.n,])
		}

#		print("init calcs...")

#		var.sites = intersect(which(parents[,1]!=parents[,2] ), which( row.sum!=0 & row.sum!=num.sam ) ) 
		## MENDELS, ERRANTS, POSITIONS by VAR SITES
		mendels.c   = mendels[coding.sites,][non.sites.c,]
		pos.seg.c = as.matrix(positions[coding.sites,1][seg.sites.c])
		pos.men.c = as.matrix(positions[coding.sites,1][non.sites.c])
		
		if (length(noncod.sites)>0)
		{
			mendels.n   = mendels[noncod.sites,][non.sites.n,]		
			pos.seg.n = as.matrix(positions[noncod.sites,1][seg.sites.n])
			pos.men.n = as.matrix(positions[noncod.sites,1][non.sites.n])
		}



#		print("more calcs...")	
		## PULL OUT CHROMO-SPECIFIC UNIQ AND COV VALUES
		chromo    = paste("MAL",j,sep ="")
		uniq.mal  = uniq[which(uniq[,1] == chromo),]
		cov.mal   = total.cov[which(total.cov[,1] == chromo),];

#		## TRANSLATE TO POSITION AND COV NUMBERING
		these.cov.seg.c = which(pos.seg.c[,1] %in% cov.mal[,2])
		those.cov.seg.c = which(cov.mal[,2] %in% pos.seg.c[,1])
		these.cov.men.c = which(pos.men.c[,1] %in% cov.mal[,2])
		those.cov.men.c = which(cov.mal[,2] %in% pos.men.c[,1])

		if (length(noncod.sites)>0)
		{
			these.cov.seg.n = which(pos.seg.n[,1] %in% cov.mal[,2])
			those.cov.seg.n = which(cov.mal[,2] %in% pos.seg.n[,1])
			these.cov.men.n = which(pos.men.n[,1] %in% cov.mal[,2])
			those.cov.men.n = which(cov.mal[,2] %in% pos.men.n[,1])
		}
		## GET THE ERR, MEN, POS and COV THAT CORRESPOND

		mendel.cov.c = mendels.c[these.cov.men.c,];
		
		cov.seg.c  = cov.mal[those.cov.seg.c,];
		cov.men.c  = cov.mal[those.cov.men.c,];
		
		mendel.cov.n = mendels.n[these.cov.men.n,];
		cov.seg.n  = cov.mal[those.cov.seg.n,];
		cov.men.n  = cov.mal[those.cov.men.n,];
	#	print("still more ...")	
		for (k in 1:length(cov.quant))
		{
			if (k ==1)
			{
				this.cov.men.c = which(cov.men.c[,3] < cov.quant[k])
				this.cov.men.n = which(cov.men.n[,3] < cov.quant[k])
			}
			else
			{
				this.cov.men.c = intersect(which(cov.men.c[,3] < cov.quant[k]), which(cov.men.c[,3] >= cov.quant[k-1]))
				this.cov.men.n = intersect(which(cov.men.n[,3] < cov.quant[k]), which(cov.men.n[,3] >= cov.quant[k-1]))
			}

			if (length(this.cov.men.c)>0)
			{
				len.men.cov.c[k] = len.men.cov.c[k] + length(this.cov.men.c);
				num.men.cov.c[k] = num.men.cov.c[k] + sum(mendel.cov.c[this.cov.men.c,]);
			}
			if ((length(noncod.sites)>0)&(length(this.cov.men.n)>0))
			{
				len.men.cov.n[k] = len.men.cov.n[k] + length(this.cov.men.n);						
				num.men.cov.n[k] = num.men.cov.n[k] + sum(mendel.cov.n[this.cov.men.n,]);
			}
		}

#		print("done");

		these.uniq.men.c = which(pos.men.c[,1] %in% uniq.mal[,2])
		those.uniq.men.c = which(uniq.mal[,2] %in% pos.men.c[,1])

		these.uniq.men.n = which(pos.men.n[,1] %in% uniq.mal[,2])
		those.uniq.men.n = which(uniq.mal[,2] %in% pos.men.n[,1])


		mendel.uniq.c = mendels.c[these.uniq.men.c,];
		unique.men.c = uniq.mal[those.uniq.men.c,];

		mendel.uniq.n = mendels.n[these.uniq.men.n,];
		unique.men.n = uniq.mal[those.uniq.men.n,];


#		print("still more ...")	
		for (k in 1:length(valid.uni))
		{
			this.uniq.men.c = which(unique.men.c[,3]==valid.uni[k])
			this.uniq.men.n = which(unique.men.n[,3]==valid.uni[k])

			if (length(this.uniq.men.c)>0)
			{
				len.men.uniq.c[k] = len.men.uniq.c[k] + length(this.uniq.men.c);
				num.men.uniq.c[k] = num.men.uniq.c[k] + sum(mendel.uniq.c[this.uniq.men.c,]);
			}
			if ((length(noncod.sites)>0)&(length(this.uniq.men.n)>0))
			{
				len.men.uniq.n[k] = len.men.uniq.n[k] + length(this.uniq.men.n);
				num.men.uniq.n[k] = num.men.uniq.n[k] + sum(mendel.uniq.n[this.uniq.men.n,]);
			}
		}
	} 

	if(r == 1)
	{
		postscript("../Figures/men.pipe.uniq.coding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,num.men.uniq.c/(20*len.men.uniq.c),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Mendelian error rate", xlab  = "Uniqueness score", ,ylim=c(0,0.08),main = "Coding")
		devs[1] = dev.cur();
		postscript("../Figures/men.pipe.uniq.noncoding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,num.men.uniq.n/(20*len.men.uniq.n),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Mendelian error rate", xlab  = "Uniqueness score", ,ylim=c(0,0.08), main = "Noncoding")
		devs[2] = dev.cur();
		postscript("../Figures/men.pipe.uniq.both.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
		plot(valid.uni,(num.men.uniq.c+num.men.uniq.n)/(20*(len.men.uniq.c+len.men.uniq.n)),type = "l",col =colors[r],pch = 19,lwd= 2,ylab = "Mendelian error rate", xlab  = "Uniqueness score", ,ylim=c(0,0.08), main = "All")
		devs[3] = dev.cur();
		postscript("../Figures/men.pipe.cov.coding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),num.men.cov.c/(20*len.men.cov.c),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.09),xlab = "Percentile coverage",ylab="Mendelian error rate",axes=FALSE,main = "Coding");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[4] = dev.cur();
		postscript("../Figures/men.pipe.cov.noncoding.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),num.men.cov.n/(20*len.men.cov.n),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.09),xlab = "Percentile coverage",ylab="Mendelian error rate",axes=FALSE, main = "Noncoding");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[5] = dev.cur();
		postscript("../Figures/men.pipe.cov.both.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);		
		plot(seq(100,1,by =-2),(num.men.cov.n+num.men.cov.c)/(20*(len.men.cov.n+len.men.cov.c)),type = "l",col =colors[r],pch = 19,lwd= 2,ylim=c(0,0.09),xlab = "Percentile coverage",ylab="Mendelian error rate",axes=FALSE, main = "All");
		axis(1,labels = seq(0,100,by = 10), at = seq(100,0,by=-10));
		axis(2,labels=TRUE);
		box()
		devs[6] = dev.cur();
		
	}
	else
	{
	#	print("adding...");
		dev.set(devs[1]);
		points(valid.uni,num.men.uniq.c/(20*len.men.uniq.c),type = "l",col =colors[r],pch = 19,lwd=2)

		if (length(noncod.sites)>0)
		{
			dev.set(devs[2]);
			points(valid.uni,num.men.uniq.n/(20*len.men.uniq.n),type = "l",col =colors[r],pch = 19,lwd=2)
		}

		dev.set(devs[3]);
		points(valid.uni,(num.men.uniq.c+num.men.uniq.n)/(20*(len.men.uniq.c+len.men.uniq.n)),type = "l",col =colors[r],pch = 19,lwd=2)

		dev.set(devs[4]);		
		points(seq(100,1,by =-2),num.men.cov.c/(20*len.men.cov.c),type = "l",col =colors[r],pch = 19,lwd=2)
	
		if (length(noncod.sites)>0)
		{
			dev.set(devs[5]);
			points(seq(100,1,by =-2),num.men.cov.n/(20*len.men.cov.n),type = "l",col =colors[r],pch = 19,lwd=2)
		}

		dev.set(devs[6]);
		points(seq(100,1,by =-2),(num.men.cov.n+num.men.cov.c)/(20*(len.men.cov.n+len.men.cov.c)),type = "l",col =colors[r],pch = 19,lwd=2)
	}
}

for (i in 1:length(devs))
{
	dev.set(devs[i]);
	legend(x ="topright", legend = cat.names, col = colors, pch = 19,cex = 1.3)
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


out.mat.n = as.matrix(rbind(num.snp.n,num.seg.n,num.non.n,num.men.n));
colnames(out.mat.n) = c("SNPs", "Seg", "Non","Men");
rownames(out.mat.n) =cat.names;
write.table(out.mat.n,"out.mat.n", quote= FALSE, sep= " & ");

for (i in 1:length(devs))
{
	dev.off(devs[i])
}
