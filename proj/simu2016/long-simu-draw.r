simu.rdata<-c("power.L1.mn.ar1", 
	"power.L2.mn.ar1",
	"power.L3.mn.ar1",
	"power.L4.mn.ar1",
	"power.L5.mn.ar1",
	"power.L6.mt.ar1",
	"power.L7.mt.ar1",
	"power.L8.msn.ar1",
	"power.L9.msn.ar1",
	"power.L10.mmn.ar1",
	"power.L11.mmn.ar1",
	"power.L12.mn.sad",
	"power.L13.mn.sad",
	"power.L14.mn.cm",
	"power.L15.mn.cm",
	"power.L16.mn.ar1");


bar.subplot<-function(p.dat, ylim=c(0,1),  main.title, test.no=0, legend=F)
{
	#if (test.no==0) 
	#	ylab<-"Power" 
	#else		
	#	ylab<-paste( "Power(G" ,test.no, ")", sep=""); 

	if (test.no==0) 
		ylab<-"" 
	else		
		ylab<-"Power";
	
	plot(1,1, type="n", xlim=c(0.4, dim(p.dat)[2]+0.6), ylim=ylim,xaxt="n",cex=0.5, cex.axis=0.75, ylab=ylab, xlab="Sample Size" )

	cols <- c("mistyrose",  "cornsilk", "lightblue", "darkblue", "gray", "black")

	p.dat <- p.dat[c(5,6,1,2,3,4), ];
	
	for(i in 1:dim(p.dat)[2])
	{
		bar.w <- 0.8/NROW(p.dat);
		for(j in 1:NROW(p.dat))
			rect(0.6+1*(i-1)+(j-1)*bar.w, 0, 0.6+1*(i-1)+j*bar.w, p.dat[j,i], col = cols[j], border = "black" );
	}
	axis(1, at=c(1:3), labels=c(500, 1000,2500), cex.axis=0.75)
	title(main.title, cex=0.7);

	if(legend)
	{
		rect( 0.5, 1.01, 0.62, 0.95, col=cols[1], border="black")
		text( 0.7-0.02, 0.98, "LSKAT", cex=0.6, adj=c(0, 0.5));

		rect( 0.5, 0.93, 0.62, 0.87, col=cols[2], border="black")
		text( 0.7-0.02, 0.90, "LBT", cex=0.6, adj=c(0, 0.5));

		rect( 0.5, 0.85, 0.62, 0.79, col=cols[3], border="black")
		text( 0.7-0.02, 0.82, "SKAT_baseline)", cex=0.6, adj=c(0, 0.5));

		rect( 0.5, 0.77, 0.62, 0.71, col=cols[4], border="black")
		text( 0.7-0.02, 0.74, "BT_baseline", cex=0.6, adj=c(0, 0.5));

		rect( 0.5, 0.69, 0.62, 0.63, col=cols[5], border="black")
		text( 0.7-0.02, 0.66, "SKAT_average", cex=0.6, adj=c(0, 0.5));

		rect( 0.5, 0.61, 0.62, 0.55, col=cols[6], border="black")
		text( 0.7-0.02, 0.58, "BT_average", cex=0.6, adj=c(0, 0.5));
	}
}

draw.bar<-function( test.no, pdf.file, loops=100)
{
	r.power0 <- c();
	if(is.numeric(test.no))
	{
		load( paste(simu.rdata[test.no], 500, "rdata", sep=".") )
		r.power0 <- r.power6;
		load( paste(simu.rdata[test.no], 1000, "rdata", sep=".") );
		r.power0 <- rbind( r.power0, r.power6) ;
		load( paste(simu.rdata[test.no], 2500, "rdata", sep=".") );
		r.power0 <- rbind( r.power0, r.power6);
	}
	else
	{
		file.rdata <- test.no;
		load( file.rdata[1] )
		r.power0 <- r.power6;
		load( file.rdata[2] );
		r.power0 <- rbind( r.power0, r.power6) ;
		load( file.rdata[3]);
		r.power0 <- rbind( r.power0, r.power6);
		
		test.no <- 1;
	}

	pp1<- r.power0[which(r.power0[,2]==1), c(6,7,8,9,3,4)];
	pp1<- t(pp1/loops);

	pp2<- r.power0[which(r.power0[,2]==2), c(6,7,8,9,3,4)];
	pp2<- t(pp2/loops);

	pp3<- r.power0[which(r.power0[,2]==3), c(6,7,8,9,3,4)];
	pp3<- t(pp3/loops);

	pdf(pdf.file, width=6, height=2.1);
	
	opar <- par(no.readonly=TRUE);
	par(pty = "m", lwd=1, oma=c(0,0,0,0) );

	#lay.mat<-rbind(
	#	c( 0,   rep(0,3),  0 ),
	#	c( 7,   1:3, 8 ),
	#	c( 7,   4:6, 8 ),
	#	c( rep(9, 5) ))
	#
	#nf <- layout(lay.mat, c(0.01, 2,  2,  2, 0.01), c(0.01, 1,  1, 0.01 ) )

	lay.mat<-rbind(
		c( 0,   rep(0,3),  0 ),
		c( 4,   1:3, 5 ),
		c( rep(6, 5) ));
	
	nf <- layout(lay.mat, c(0.01, 2,  2,  2, 0.01), c(0.01, 1, 0.01 ) );

	par(mar=c( 3, 3, 2.5, 0.5 ));

	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(pp3, ylim=c(0,1), "50% positive", test.no, legend=T)

	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(pp2, ylim=c(0,1), "80% positive");

	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(pp1, ylim=c(0,1), "100% positive");

	dev.off();
}

draw.bar( 1, "Power-L1.pdf")
draw.bar( 2, "Power-L2.pdf")
draw.bar( 3, "Power-L3.pdf")
draw.bar( 4, "Power-L4.pdf")
draw.bar( 5, "Power-L5.pdf")
draw.bar( 6, "Power-L6.pdf")
draw.bar( 7, "Power-L7.pdf")
draw.bar( 8, "Power-L8.pdf")
draw.bar( 9, "Power-L9.pdf")
draw.bar( 10, "Power-L10.pdf")
draw.bar( 11, "Power-L11.pdf")
draw.bar( 12, "Power-L12.pdf")
draw.bar( 13, "Power-L13.pdf")
draw.bar( 14, "Power-L14.pdf")
draw.bar( 15, "Power-L15.pdf")
draw.bar( 16, "Power-L16.pdf")



draw.bar( c("power.L1.mn.ar1.500.t2.rdata", "power.L1.mn.ar1.1000.t2.rdata", "power.L1.mn.ar1.2500.t2.rdata"), "Power-L1-t2.pdf")
