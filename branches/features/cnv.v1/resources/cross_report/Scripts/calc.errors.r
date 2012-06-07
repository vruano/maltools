cat("-- reading arguments\n", sep = "");
cmd_args = commandArgs();

prefix = cmd_args[4];
figure = cmd_args[5];
cross_1 = cmd_args[6];
cross_2 = cmd_args[7];

prefix.parents = paste(prefix,"parents.MAL",sep="");
prefix.progeny = paste(prefix,"progeny.MAL",sep="");
prefix.paths = paste(prefix,"paths.MAL",sep="");
prefix.positions = paste(prefix,"positions.MAL",sep="");

for (i in 1:14)
{
	mendel = paste(prefix,"mendels.MAL",sep="");
	mendel = paste(mendel,i,sep="");
	mendel = paste(mendel,".txt",sep="");

	errant = paste(prefix,"errant.MAL",sep="");
	errant = paste(errant,i,sep="");
	errant = paste(errant,".txt",sep="");

	file.name = paste(figure,"cross",sep = "");
	file.name = paste(file.name,".MAL",sep="");
	file.name = paste(file.name,i,sep="");
	file.name = paste(file.name,".ps",sep="");

	parent.file = paste(prefix.parents,i,sep="");

        print(paste("Readding parents... ",parent.file));	
	parents = read.table(parent.file,header = FALSE,sep =",");

	progeny.file = paste(prefix.progeny,i,sep="");
	progeny = read.table(progeny.file, header = FALSE, sep = ",");

	num.samples = dim(progeny)[2];

	paths.file = paste(prefix.paths,i,sep="");
	print(paths.file)
	paths = read.table(paths.file,header = FALSE,sep=",");


	positions.file = paste(prefix.positions,i,sep="");
	print(positions.file)
	positions = read.table(positions.file,header = FALSE,sep=",");
	positions = positions[,1];

	postscript(file.name,paper="special", width = 14, height = 10,horizontal = FALSE)
	layout(matrix(c(1,1,2,2),byrow = TRUE, ncol = 2),heights = c(3,1));
	par(mar=c(0,7,1.5,2),oma=c(0,0,2,0))
	out.mat = cbind(t(paths),(parents))
	image(as.matrix(out.mat),col = grey((0:10)/10),axes= FALSE);

	mendel.errors = matrix(rep(0,num.samples*length(positions)),ncol = length(positions));


	joint.snps = which(parents[,2] != parents[,1]);

	for (k in 1:num.samples)
	{
		mendel.errors[k,which((progeny[,k] != parents[,1])&(progeny[,k] != parents[,2]))] = 1;
	}


	for (k in 1:num.samples)
	{
		num.mendel = sum(mendel.errors[k,])
		if (num.mendel > 0)
		{
			x = which(mendel.errors[k,]==1)/length(mendel.errors[k,]);
			y = rep(0.002+(k-1)/(num.samples+1),num.mendel);
	
			points( x, y,col ="red",pch = 19, cex=0.5)
		}
	}

	for(k in 0:(num.samples))
	{
		if (i  == num.samples-1)
		{
			lines(x= c(0,1), y = c(0.024+k/(num.samples+1),0.024+k/(num.samples+ 1)),col = "purple",lwd =6,cex = 1)
		}
		else
		{
			lines(x= c(0,1), y = c(0.024+k/(num.samples+1),0.024+k/(num.samples+1)),col = grey(0.5),lwd =1,cex = 1)
		}
	}

	num.positions = length(positions)
	hmm.errors = matrix(rep(0,num.samples*length(positions)),ncol = length(positions));

	for (k in 1:num.samples)
	{
		p.1 = rep(0,num.positions);
		noughts = which(paths[k,]==0);
		calls.hmm = union(which(parents[noughts,1]!= progeny[noughts,k]),which(parents[-noughts,2]!=progeny[-noughts,k]))
	  	hmm.errors[k,calls.hmm] = 1;
		hmm.errors[k,-joint.snps]=0;

		if (length(calls.hmm)>0)
		{
	 		points( x = calls.hmm/length(calls.hmm), y = rep(0.009+(k-1)/(num.samples+1),length(calls.hmm)),col ="blue",pch = 19, cex =0.5)
		}	
	}

        print(paste("Wrinting errants ",errant))
	write.table(t(hmm.errors), file=errant, col.names = FALSE, row.names = FALSE, sep = "\t")
        print(paste("Wrinting mendel ",mendel))
	write.table(t(mendel.errors), file=mendel, col.names = FALSE, row.names = FALSE, sep = "\t")

	total.num = num.samples + 2;
	axis(2, labels = c(seq(1,num.samples,by = 1),cross_1,cross_2), at = seq(0.5,total.num+0.5,length.out= (total.num))/total.num, cex.axis = 0.8);
	box()
	string = paste("HMM path, with Mendelian errors and errant SNP calls: Chromosome ",i,sep = "")
	title(main = string, lty = 3)

	snp.positions = positions;
	num.snps = length(snp.positions);

	max.val = max(snp.positions);
	min.val = min(snp.positions);
	plot(seq(0,1,length.out = max.val),rep(1,max.val), type="n", axes=F, xlab= "SNP Position", ylab = "",xlim = c(0.035,0.96), ylim = c(0.03,1))
	par(mar=c(4,7,2,2),oma=c(1,0,0,0))
	box();

	for (k in seq(1,num.snps,by = 20))
	{
		lines(x = c((snp.positions[k])/max.val,(k-1)/num.snps),y=c(0,1),col ="black");
	}

	axis(1,at = c(0.04,0.25,0.5,0.75,0.96), label = round(c(0.04,0.25,0.5,0.75,0.96)*max.val))
	title( xlab = "SNP position")

	dev.off();

}





