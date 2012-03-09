#prefixes = c("../GATK/3D7.HB3","../BWA2SOM/3D7.HB3","Data.Type/3D7.HB3","../CORTEX/3D7.HB3");
#cat.names =c("GATK","BWA2SOM","PGV","CORTEX");
cat.names = c("Original", "OK", "Coverage", "WCSA", "Biallelic", "Uniqueness", "Typable")
prefixes = c("Data.Orig/3D7.HB3","Data.OK/3D7.HB3","Data.Cov/3D7.HB3","Data.WCSA/3D7.HB3","Data.Bi/3D7.HB3","Data.Uniq/3D7.HB3","Data.Type/3D7.HB3");
uniq = read.table("uniq.data",header = FALSE,sep = "\t");
cov   = read.table("original.coverage",header =TRUE, sep ="\t");
total.cov = rowSums(cov[,4:7]);
total.cov = cbind(cov[,1:2],total.cov);
#rm("cov");

cov.quant = quantile(total.cov[,3],seq(0.05,1,by =.05));
num.cat= length(prefixes);
colors = palette(rainbow(num.cat, s = 0.5, v = 0.75))

cov.seq = seq(0.02,1,by=0.02)
cov.quant = quantile(total.cov[,3],cov.seq,names = FALSE);


num.snp  = rep(0,num.cat);
num.seg  = rep(0,num.cat);
num.non  = rep(0,num.cat);
num.men  = rep(0,num.cat);
num.err  = rep(0,num.cat);

postscript("../Figures/test.pipes.coding.uniq.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/men.sites.uniq.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/men.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.men.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.men.sites.uniq.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);
#postscript("../Figures/err.sites.cov.pgv.ps",horizontal = FALSE,paper = "special", height = 7, width = 10);

cnc = rep(list(),14)
for( r in 1)
{
	cnc[[i]] = read.table(paste("CNC/cnc",i,sep=""),header = FALSE, sep="");
}

for (r in 1)
{
	prefix = prefixes[r];

	print(r);
	num.men.cov = rep(0,length(cov.seq));
	num.err.cov = rep(0,length(cov.seq));
	num.tot.cov = rep(0,length(cov.seq));
	len.seg.cov = rep(0,length(cov.seq));
	len.men.cov = rep(0,length(cov.seq));

	num.men.uniq = rep(0,14);
	num.err.uniq = rep(0,14);
	num.tot.uniq = rep(0,14);
	len.seg.uniq = rep(0,14);
	len.men.uniq = rep(0,14);

	for (j in 1)
	{

		## FILE NAMES
		if ((prefix == "../GATK/3D7.HB3")|(prefix == "../BWA2SOM/3D7.HB3"))
		{
			errant.file = paste(prefix,"errant.MAL.",sep="/");
			errant.file = paste(errant.file,j,sep="");
			errant.file = paste(errant.file,".txt",sep="");

			mendel.file = paste(prefix,"mendels.MAL.",sep="/");
			mendel.file = paste(mendel.file,j,sep="");
			mendel.file = paste(mendel.file,".txt",sep="");
		}
		else
		{
			errant.file = paste(prefix,"errant.MAL",sep="/");
			errant.file = paste(errant.file,j,sep="");
			errant.file = paste(errant.file,".txt",sep="");

			mendel.file = paste(prefix,"mendels.MAL",sep="/");
			mendel.file = paste(mendel.file,j,sep="");
			mendel.file = paste(mendel.file,".txt",sep="");

		}

		positions.file = paste(prefix,"positions.MAL",sep="/");
		positions.file = paste(positions.file,j,sep="");

		parents.file = paste(prefix,"parents.MAL",sep="/");
		parents.file = paste(parents.file,j,sep="");

		progeny.file = paste(prefix,"progeny.MAL",sep="/");
		progeny.file = paste(progeny.file,j,sep="");

		mendels   = read.table(mendel.file,header =FALSE, sep ="\t");
		errants   = read.table(errant.file,header =FALSE, sep ="\t");

		positions = read.table(positions.file);

		parents   = read.table(parents.file,header = FALSE, sep =",");
		progeny   = read.table(progeny.file,header = FALSE, sep =",");


		## TOTAL DATA, NUM SAMPLES, ROW SUMS
		totals    = cbind(parents,progeny);
		num.sam   = dim(totals)[2];
		row.sum   = rowSums(totals); 
		
		## FIND VARIABLE SITES 
		seg.sites = which(parents[,1] != parents[,2]);
		non.sites = which(parents[,1] == parents[,2]);
#		non.sites = intersect(non.sites, which( row.sum!=0 & row.sum!=num.sam ) ) 

		## CALCULATE 
		num.snp[r]  = num.snp[r]  + dim(parents)[1];
		print(num.snp[r])
		num.seg[r]  = num.seg[r]  + length(seg.sites);
		num.non[r]  = num.non[r]  + length(non.sites);
		num.men[r]  = num.men[r]  + sum(mendels[non.sites,])
		num.err[r]  = num.err[r]  + sum(errants[non.sites,])
	

#		var.sites = intersect(which(parents[,1]!=parents[,2] ), which( row.sum!=0 & row.sum!=num.sam ) ) 
		## MENDELS, ERRANTS, POSITIONS by VAR SITES
		mendels   = mendels[non.sites,]
		errants   = errants[non.sites,]
		pos.seg = as.matrix(positions[seg.sites,])
		pos.men = as.matrix(positions[non.sites,])

		#print(c(length(pos.seg),length(pos.men)));
		##  GET COUNTS FOR 
	
		## PULL OUT CHROMO-SPECIFIC UNIQ AND COV VALUES
		chromo    = paste("MAL",j,sep ="")
		uniq.mal  = uniq[which(uniq[,1] == chromo),]
		cov.mal   = total.cov[which(total.cov[,1] == chromo),];

		## TRANSLATE TO POSITION AND COV NUMBERING
		these.cov.seg = which(pos.seg[,1] %in% cov.mal[,2])
		those.cov.seg = which(cov.mal[,2] %in% pos.seg[,1])
		these.cov.men = which(pos.men[,1] %in% cov.mal[,2])
		those.cov.men = which(cov.mal[,2] %in% pos.men[,1])


		## GET THE ERR, MEN, POS and COV THAT CORRESPOND
		errant.cov = errants[these.cov.men,];
		mendel.cov = mendels[these.cov.men,];

		cov.seg    = cov.mal[those.cov.seg,];
		cov.men    = cov.mal[those.cov.men,];

		for (k in 1:length(cov.quant))
		{

			if (k ==1)
			{
				this.cov.seg = which(cov.seg[,3] < cov.quant[k])
				this.cov.men = which(cov.men[,3] < cov.quant[k])
			}
			else
			{
				this.cov.seg = intersect(which(cov.seg[,3] < cov.quant[k]), which(cov.seg[,3] >= cov.quant[k-1]))
				this.cov.men = intersect(which(cov.men[,3] < cov.quant[k]), which(cov.men[,3] >= cov.quant[k-1]))
			}
			#print(length(this.cov.seg));
			if (length(this.cov.seg)>0)
			{
				len.seg.cov[k] = len.seg.cov[k] + length(this.cov.seg);
				num.err.cov[k] = num.err.cov[k] + sum(errant.cov[this.cov.men,]);				
			}
			if (length(this.cov.men)>0)
			{
				len.men.cov[k] = len.men.cov[k] + length(this.cov.men);			
				num.men.cov[k] = num.men.cov[k] + sum(mendel.cov[this.cov.men,]);
			}

			num.tot.cov[k] = num.tot.cov[k] + num.men.cov[k] + num.err.cov[k];
		}
		#print("done");
		these.uniq.seg = which(pos.seg[,1] %in% uniq.mal[,2])
		those.uniq.seg = which(uniq.mal[,2] %in% pos.seg[,1])
		these.uniq.men = which(pos.men[,1] %in% uniq.mal[,2])
		those.uniq.men = which(uniq.mal[,2] %in% pos.men[,1])


		errant.uniq = errants[these.uniq.men,];
		mendel.uniq = mendels[these.uniq.men,];

		unique.seg = uniq.mal[those.uniq.seg,];
		unique.men = uniq.mal[those.uniq.men,];

		valid.uni = seq(18,44,by=2);

		for (k in 1:length(valid.uni))
		{

			this.uniq.seg = which(unique.seg[,3]==valid.uni[k])
			this.uniq.men = which(unique.men[,3]==valid.uni[k])

			if (length(this.uniq.seg)>0)
			{
				len.seg.uniq[k] = len.seg.uniq[k] + length(this.uniq.seg);
				num.err.uniq[k] = num.err.uniq[k] + sum(errant.uniq[this.uniq.men,]);	
			}
			if (length(this.uniq.men)>0)
			{
				len.men.uniq[k] = len.men.uniq[k] + length(this.uniq.men);
				num.men.uniq[k] = num.men.uniq[k] + sum(mendel.uniq[this.uniq.men,]);
			}
			num.tot.uniq[k] = num.tot.uniq[k] + num.men.uniq[k] + num.err.uniq[k];
		}
	} 

	if(r == 1)
	{
		plot(valid.uni,num.err.uniq/(20*len.men.uniq),type = "l",col =colors[r],pch = 19,lwd= 3,ylab = "Errant call rate", xlab  = "Uniqueness score", ,ylim=c(0,0.033))
#		plot((cov.seq)*100,num.err.cov[length(cov.seq):1]/(20*len.men.cov[length(cov.seq):1]),type = "l",col =colors[r],pch = 19,lwd= 3,ylim=c(0,0.019),xlab = "Percentile of coverage",ylab="Total error rate")
	}
	else
	{
		points(valid.uni,num.err.uniq/(20*len.men.uniq),type = "l",col =colors[r],pch = 19,lwd= 3,ylim=c(0,0.15))
#		points((cov.seq)*100,num.err.cov[length(cov.seq):1]/(20*len.men.cov[length(cov.seq):1]),type = "l",col =colors[r],pch = 19, lwd= 3)
	}
	legend(x ="topleft", legend =cat.names, col = colors, pch = 19,cex = 1.3)
#	legend(x ="topleft", legend = c("Original","OK","Good coverage","WCSA","Biallelic","Uniqueness", "Typable"), col = colors, pch = 19,cex = 1.3)
}
#out.mat = as.matrix(rbind(num.snp,num.seg,num.non,num.err,num.men));
##rownames(out.mat) = c("SNPs", "Seg", "Non","Err","Men");
#colnames(out.mat) =cat.names;
dev.off()
