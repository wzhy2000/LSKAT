if(0)
{
simu.rdata<-c("L1.mn.ar1",
	"L2.mn.ar1",
	"L3.mn.ar1",
	"L4.mn.ar1",
	"L5.mn.ar1",
	"L6.mt.ar1",
	"L7.mt.ar1",
	"L8.msn.ar1",
	"L9.msn.ar1",
	"L10.mmn.ar1",
	"L11.mmn.ar1",
	"L12.mn.sad",
	"L13.mn.sad",
	"L14.mn.cm",
	"L15.mn.cm",
	"L16.mn.ar1");

simu.folder<-c("L1-mn-ar1",
	"L2-mn-ar1",
	"L3-mn-ar1",
	"L4-mn-ar1",
	"L5-mn-ar1",
	"L6-mt-ar1",
	"L7-mt-ar1",
	"L8-msn-ar1",
	"L9-msn-ar1",
	"L10-mmn-ar1",
	"L11-mmn-ar1",
	"L12-mn-sad",
	"L13-mn-sad",
	"L14-mn-cm",
	"L15-mn-cm",
	"L16-mn-ar1");
}


simu.rdata<-c("L1.1k.mn.ar1",
	"L4.1k.mn.ar1",
	"L5.1k.mn.ar1",
	"L12.1k.mn.sad",
	"L14.1k.mn.cm");

simu.folder<-c("L1-mn-ar1",
	"L4-mn-ar1",
	"L5-mn-ar1",
	"L12-mn-sad",
	"L14-mn-cm");


check_type1_ret<-function( simu.folder, simu.rdata)
{
	unfind <-c();

	r500 <-c();
	for(i in 1:100)
	{
		rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".500.", i ,".rdata", sep="");
		x<- try( load( rdata.file), TRUE ) ;
		if(class(x)!="try-error")
			r500 <- rbind( r500, r.1 )
		else
		{
			cat("NOT Find ", rdata.file, "\n");
			unfind<-c(unfind, rdata.file);
		}
	}


	r1000 <-c();
	for(i in 1:100)
	{
		rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".1000.", i ,".rdata", sep="");
		x<- try( load( rdata.file), TRUE ) ;
		if(class(x)!="try-error")
			r1000 <- rbind( r1000, r.1 )
		else
		{
			cat("NOT Find ", rdata.file, "\n");
			unfind<-c(unfind, rdata.file);
		}
	}

	r2500 <-c();
	for(i in 1:100)
	{
		rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".2500.", i ,".rdata", sep="");
		x<- try( load( rdata.file), TRUE ) ;
		if(class(x)!="try-error")
			r2500 <- rbind( r2500, r.1 )
		else
		{
			cat("NOT Find ", rdata.file, "\n");
			unfind<-c(unfind, rdata.file);
		}
	}


	a5.level<-10^(-5);
	a4.level<-10^(-4);
	a3.level<-10^(-3);
	a2.level<-10^(-2);

	nsample <- 1000
	nloop <- 1000

	r.1<- r500;

	r.type1.500<-c(nloop, nsample, 500,
			  length(which( r.1[,3]<a2.level)),
			  length(which( r.1[,3]<a3.level)),
			  length(which( r.1[,3]<a4.level)),
			  length(which( r.1[,3]<a5.level)),
			  mean(r.1[,6]), mean(r.1[,7]) );


	r.1<- r1000;

	r.type1.1000<-c(nloop, nsample, 1000,
			  length(which( r.1[,3]<a2.level)),
			  length(which( r.1[,3]<a3.level)),
			  length(which( r.1[,3]<a4.level)),
			  length(which( r.1[,3]<a5.level)),
			  mean(r.1[,6]), mean(r.1[,7]) );

	r.1<- r2500;

	r.type1.2500<-c(nloop, nsample, 2500,
			  length(which( r.1[,3]<a2.level)),
			  length(which( r.1[,3]<a3.level)),
			  length(which( r.1[,3]<a4.level)),
			  length(which( r.1[,3]<a5.level)),
			  mean(r.1[,6]), mean(r.1[,7]) );

	r.type1 <- rbind(r.type1.500, r.type1.1000, r.type1.2500);

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata, ".all.rdata", sep="");
	save(par, r500, r1000, r2500, r.type1, file=rdata.file )

	show(r.type1);

	return(list(r.type1=r.type1, unfind=unfind));

}

a.level<- 10^(-6);

check_power0_ret<-function( simu.folder, simu.rdata)
{
	unfind <- c();
	r.power0 <- c();
	snp1 <- c();
	rare1 <- c();

	p.1.500<-c();
	p.2.500<-c();
	p.3.500<-c();
	p.1.1000<-c();
	p.2.1000<-c();
	p.3.1000<-c();
	p.1.2500<-c();
	p.2.2500<-c();
	p.3.2500<-c();

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".500.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.500 <- r.1; snp1<-c(snp1, mean(p.1.500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.1.500[,"lskat.snp.rare"]));
		p.2.500 <- r.2; snp1<-c(snp1, mean(p.2.500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.2.500[,"lskat.snp.rare"]));
		p.3.500 <- r.3; snp1<-c(snp1, mean(p.3.500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.3.500[,"lskat.snp.rare"]));

		r.power <- c();
		r.power<-rbind( r.power, c(500, 100, length(which(r.1[,13]<a.level)),
										     length(which(r.1[,14]<a.level)),
										     length(which(r.1[,18]<a.level)),
										     length(which(r.1[,19]<a.level)),
										     length(which(r.1[,20]<a.level)),
										     length(which(r.1[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 80,  length(which(r.2[,13]<a.level)),
										     length(which(r.2[,14]<a.level)),
										     length(which(r.2[,18]<a.level)),
										     length(which(r.2[,19]<a.level)),
										     length(which(r.2[,20]<a.level)),
										     length(which(r.2[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 50,  length(which(r.3[,13]<a.level)),
										     length(which(r.3[,14]<a.level)),
										     length(which(r.3[,18]<a.level)),
										     length(which(r.3[,19]<a.level)),
										     length(which(r.3[,20]<a.level)),
										     length(which(r.3[,21]<a.level)) ) );

		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".1000.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.1000 <- r.1; snp1<-c(snp1, mean(p.1.1000[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.1.1000[,"lskat.snp.rare"]));
		p.2.1000 <- r.2; snp1<-c(snp1, mean(p.2.1000[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.2.1000[,"lskat.snp.rare"]));
		p.3.1000 <- r.3; snp1<-c(snp1, mean(p.3.1000[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.3.1000[,"lskat.snp.rare"]));

		r.power <- c();

		r.power<-rbind( r.power, c(500*2, 100, length(which(r.1[,13]<a.level)),
										     length(which(r.1[,14]<a.level)),
										     length(which(r.1[,18]<a.level)),
										     length(which(r.1[,19]<a.level)),
										     length(which(r.1[,20]<a.level)),
										     length(which(r.1[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500*2, 80,  length(which(r.2[,13]<a.level)),
										     length(which(r.2[,14]<a.level)),
										     length(which(r.2[,18]<a.level)),
										     length(which(r.2[,19]<a.level)),
										     length(which(r.2[,20]<a.level)),
										     length(which(r.2[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500*2, 50,  length(which(r.3[,13]<a.level)),
										     length(which(r.3[,14]<a.level)),
										     length(which(r.3[,18]<a.level)),
										     length(which(r.3[,19]<a.level)),
										     length(which(r.3[,20]<a.level)),
										     length(which(r.3[,21]<a.level)) ) );


		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".2500.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.2500 <- r.1; snp1<-c(snp1, mean(p.1.2500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.1.2500[,"lskat.snp.rare"]));
		p.2.2500 <- r.2; snp1<-c(snp1, mean(p.2.2500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.2.2500[,"lskat.snp.rare"]));
		p.3.2500 <- r.3; snp1<-c(snp1, mean(p.3.2500[,"lskat.snp.total"])); rare1<-c(rare1, mean(p.3.2500[,"lskat.snp.rare"]));

		r.power<-rbind( r.power, c(500*5, 100, length(which(r.1[,13]<a.level)),
										     length(which(r.1[,14]<a.level)),
										     length(which(r.1[,18]<a.level)),
										     length(which(r.1[,19]<a.level)),
										     length(which(r.1[,20]<a.level)),
										     length(which(r.1[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500*5, 80,  length(which(r.2[,13]<a.level)),
										     length(which(r.2[,14]<a.level)),
										     length(which(r.2[,18]<a.level)),
										     length(which(r.2[,19]<a.level)),
										     length(which(r.2[,20]<a.level)),
										     length(which(r.2[,21]<a.level)) ) );
		r.power<-rbind( r.power, c(500*5, 50,  length(which(r.3[,13]<a.level)),
										     length(which(r.3[,14]<a.level)),
										     length(which(r.3[,18]<a.level)),
										     length(which(r.3[,19]<a.level)),
										     length(which(r.3[,20]<a.level)),
										     length(which(r.3[,21]<a.level)) ) );

		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	r.power0<-cbind(r.power0, SNPs=snp1, rare=rare1);

	r.power0 <- r.power0[,c(1,2,9,10,3,4,5,6,7,8 )];
	colnames(r.power0) <-c("Sample", "Model", "Total SNP", "Rare SNP", "LSKAT", "LBL", "SKAT-BL", "BT-BL", "SKAT-MU", "BT-MU" );

	rdata.file <- paste(simu.folder, "/power0.", simu.rdata, ".all.rdata", sep="");
	save(par, r.power0, p.1.500, p.2.500, p.3.500, p.1.1000, p.2.1000, p.3.1000, p.1.2500, p.2.2500, p.3.2500, file=rdata.file )

	show( r.power0 );

	return(list(r.power0=r.power0, unfind=unfind));
}

check_power2_ret<-function( simu.folder, simu.rdata)
{
	unfind <- c();
	r.power0 <- c();
	snp1 <- c();
	rare1 <- c();

	p.1.500<-c();
	p.2.500<-c();
	p.3.500<-c();
	p.4.500<-c();
	p.5.500<-c();
	p.6.500<-c();

	p.1.1000<-c();
	p.2.1000<-c();
	p.3.1000<-c();
	p.4.1000<-c();
	p.5.1000<-c();
	p.6.1000<-c();

	p.1.2500<-c();
	p.2.2500<-c();
	p.3.2500<-c();
	p.4.2500<-c();
	p.5.2500<-c();
	p.6.2500<-c();

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".500.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.500 <- r.1; snp1<-c(snp1, mean(p.1.500[,8])); rare1<-c(rare1, mean(p.1.500[,9]));
		p.2.500 <- r.2; snp1<-c(snp1, mean(p.2.500[,8])); rare1<-c(rare1, mean(p.2.500[,9]));
		p.3.500 <- r.3; snp1<-c(snp1, mean(p.3.500[,8])); rare1<-c(rare1, mean(p.3.500[,9]));
		p.4.500 <- r.4; snp1<-c(snp1, mean(p.4.500[,8])); rare1<-c(rare1, mean(p.4.500[,9]));
		p.5.500 <- r.5; snp1<-c(snp1, mean(p.5.500[,8])); rare1<-c(rare1, mean(p.5.500[,9]));
		p.6.500 <- r.6; snp1<-c(snp1, mean(p.6.500[,8])); rare1<-c(rare1, mean(p.6.500[,9]));

		r.power <- c();
		r.power<-rbind( r.power, c(500, 1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)), length(which(r.1[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)), length(which(r.2[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)), length(which(r.3[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)), length(which(r.4[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)), length(which(r.5[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(500, 6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)), length(which(r.6[,5]<a.level)) ) );

		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".1000.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.1000 <- r.1; snp1<-c(snp1, mean(p.1.1000[,8])); rare1<-c(rare1, mean(p.1.1000[,9]));
		p.2.1000 <- r.2; snp1<-c(snp1, mean(p.2.1000[,8])); rare1<-c(rare1, mean(p.2.1000[,9]));
		p.3.1000 <- r.3; snp1<-c(snp1, mean(p.3.1000[,8])); rare1<-c(rare1, mean(p.3.1000[,9]));
		p.4.1000 <- r.4; snp1<-c(snp1, mean(p.4.1000[,8])); rare1<-c(rare1, mean(p.4.1000[,9]));
		p.5.1000 <- r.5; snp1<-c(snp1, mean(p.5.1000[,8])); rare1<-c(rare1, mean(p.5.1000[,9]));
		p.6.1000 <- r.6; snp1<-c(snp1, mean(p.6.1000[,8])); rare1<-c(rare1, mean(p.6.1000[,9]));

		r.power <- c();
		r.power<-rbind( r.power, c(1000,  1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)), length(which(r.1[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(1000,  2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)), length(which(r.2[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(1000,  3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)), length(which(r.3[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(1000,  4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)), length(which(r.4[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(1000,  5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)), length(which(r.5[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(1000,  6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)), length(which(r.6[,5]<a.level)) ) );

		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	rdata.file <- paste(simu.folder, "/power.", simu.rdata,".2500.rdata", sep="");
	x<-try( load(rdata.file), TRUE )
	if (class(x)!="try-error")
	{
		p.1.2500 <- r.1; snp1<-c(snp1, mean(p.1.2500[,8])); rare1<-c(rare1, mean(p.1.2500[,9]));
		p.2.2500 <- r.2; snp1<-c(snp1, mean(p.2.2500[,8])); rare1<-c(rare1, mean(p.2.2500[,9]));
		p.3.2500 <- r.3; snp1<-c(snp1, mean(p.3.2500[,8])); rare1<-c(rare1, mean(p.3.2500[,9]));
		p.4.2500 <- r.4; snp1<-c(snp1, mean(p.4.2500[,8])); rare1<-c(rare1, mean(p.4.2500[,9]));
		p.5.2500 <- r.5; snp1<-c(snp1, mean(p.5.2500[,8])); rare1<-c(rare1, mean(p.5.2500[,9]));
		p.6.2500 <- r.6; snp1<-c(snp1, mean(p.6.2500[,8])); rare1<-c(rare1, mean(p.6.2500[,9]));

		r.power <- c();
		r.power<-rbind( r.power, c(2500, 1, length(which(r.1[,2]<a.level)), length(which(r.1[,3]<a.level)), length(which(r.1[,4]<a.level)), length(which(r.1[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(2500, 2, length(which(r.2[,2]<a.level)), length(which(r.2[,3]<a.level)), length(which(r.2[,4]<a.level)), length(which(r.2[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(2500, 3, length(which(r.3[,2]<a.level)), length(which(r.3[,3]<a.level)), length(which(r.3[,4]<a.level)), length(which(r.3[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(2500, 4, length(which(r.4[,2]<a.level)), length(which(r.4[,3]<a.level)), length(which(r.4[,4]<a.level)), length(which(r.4[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(2500, 5, length(which(r.5[,2]<a.level)), length(which(r.5[,3]<a.level)), length(which(r.5[,4]<a.level)), length(which(r.5[,5]<a.level)) ) );
		r.power<-rbind( r.power, c(2500, 6, length(which(r.6[,2]<a.level)), length(which(r.6[,3]<a.level)), length(which(r.6[,4]<a.level)), length(which(r.6[,5]<a.level)) ) );

		r.power0 <- rbind(r.power0, r.power);
	}
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	r.power0<-cbind(r.power0, SNPs=snp1, rare=rare1);

	rdata.file <- paste(simu.folder, "/power0.", simu.rdata, ".all.rdata", sep="");
	save(par, r.power0, p.1.500, p.2.500, p.3.500, p.1.1000, p.2.1000, p.3.1000, p.1.2500, p.2.2500, p.3.2500, file=rdata.file )

	show( r.power0 );

	return(list(r.power0=r.power0, unfind=unfind));
}


check_type1_ret_new<-function( simu.folder, simu.rdata)
{
	unfind <-c();
	r.type <- c();

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".500.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)!="try-error")
		r.type <- rbind( r.type, r.type1)
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}


	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".1000.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)!="try-error")
		r.type <- rbind( r.type, r.type1)
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".2500.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)!="try-error")
		r.type <- rbind( r.type, r.type1)
	else
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}

	show(r.type);

	return(list(r.type1=r.type, unfind=unfind));

}

check_type1_ret_new<-function( simu.folder, simu.rdata, nsample= 100, nloop = 100)
{
	cat("***", simu.folder, simu.rdata, "\n");
	a5.level<-10^(-5);
	a4.level<-10^(-4);
	a3.level<-10^(-3);
	a2.level<-10^(-2);

	unfind <-c();
	r500 <- r1000 <- r2500 <- c();

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".500.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)=="try-error")
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}
	else
	{
		r500 <- do.call("rbind", lapply(1:2, function(i){
				  rx <- c(nloop, nsample, 500, i, mean(r.1[,7]), mean(r.1[,8]),
				  length(which( r.1[,4+i]<a2.level)),
				  length(which( r.1[,4+i]<a3.level)),
				  length(which( r.1[,4+i]<a4.level)),
				  length(which( r.1[,4+i]<a5.level)) );
			return(rx); } ) );

		r.1 <- NULL;
	}

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".1000.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)=="try-error")
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}
	else
	{
		r1000 <- do.call("rbind", lapply(1:2, function(i){
			rx <- c(nloop, nsample, 1000, i, mean(r.1[,7]), mean(r.1[,8]),
				  length(which( r.1[,4+i]<a2.level)),
				  length(which( r.1[,4+i]<a3.level)),
				  length(which( r.1[,4+i]<a4.level)),
				  length(which( r.1[,4+i]<a5.level)) );
			return(rx); } ) );
		r.1 <- NULL;
	}

	rdata.file <- paste(simu.folder, "/type1.", simu.rdata,".2500.rdata", sep="");
	x<- try( load( rdata.file), TRUE ) ;
	if(class(x)=="try-error")
	{
		cat("NOT Find ", rdata.file, "\n");
		unfind<-c(unfind, rdata.file);
	}
	else
	{
		r2500 <- do.call("rbind", lapply(1:2, function(i){
			rx <- c(nloop, nsample, 2500, i, mean(r.1[,7]), mean(r.1[,8]),
				  length(which( r.1[,4+i]<a2.level)),
				  length(which( r.1[,4+i]<a3.level)),
				  length(which( r.1[,4+i]<a4.level)),
				  length(which( r.1[,4+i]<a5.level)) );
			return(rx); } ));
		r.1 <- NULL;
	}

	r.type1 <- rbind(r500, r1000, r2500);
	r.type1[,c(7:10)] <- r.type1[,c(7:10)] /nsample/nloop;
	colnames(r.type1) <-c("Loop", "repi", "sample", "LSKAT1_LBT2", "total", "rare", "a2", "a3", "a4", "a5");
	rdata.file <- paste(simu.folder, "/type1.", simu.rdata, ".all.rdata", sep="");
	save(par, r500, r1000, r2500, r.type1, file=rdata.file )


	show(r.type1);

	return(list(r.type1=r.type1, unfind=unfind));
}


for(i in 1:16)
{
	#cat("Test Plan = ", simu.folder[i], "=======================================>\n");
	#check_type1_ret( simu.folder[i], simu.rdata[i]);
	#check_power0_ret(simu.folder[i], simu.rdata[i]);
}

for(i in 1:16)
{
	#cat("Test Plan = ", simu.folder[i], "=======================================>\n");
	#check_power2_ret(simu.folder[i], simu.rdata[i]);
}

for(i in 1:16)
{
	#cat("Test Plan = ", simu.folder[i], "=======================================>\n");
	#check_power0_ret(".", simu.rdata[i]);
}


for(i in 1:16)
{
	#cat("Test Plan = ", simu.folder[i], "=======================================>\n");
	check_type1_ret_new(".", simu.rdata[i], nsample= 1000, nloop = 1000);
}

