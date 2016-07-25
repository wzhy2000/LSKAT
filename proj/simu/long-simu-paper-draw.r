simu.rdata<-c("power.L1.1k.mn.ar1", 
	"power.L2.mn.ar1",
	"power.L3.mn.ar1",
	"power.L4.1k.mn.ar1",
	"power.L5.1k.mn.ar1",
	"power.L6.mt.ar1",
	"power.L7.mt.ar1",
	"power.L8.msn.ar1",
	"power.L9.msn.ar1",
	"power.L10.mmn.ar1",
	"power.L11.mmn.ar1",
	"power.L12.1k.mn.sad",
	"power.L13.mn.sad",
	"power.L14.1k.mn.cm",
	"power.L15.mn.cm",
	"power.L16.mn.ar1");


bar.subplot<-function(p.dat, ylim=c(0, 1.1), letter, test.no=0, legend=F,xlab=F, ylab=F)
{
	plot(1,1, type="n", xlim=c(0.4, dim(p.dat)[2]+0.6), ylim=ylim,xaxt="n",cex=0.5, cex.axis=0.75, ylab=ifelse(ylab, "Power", ""),  xlab=ifelse(xlab, "Sample Size", "") )

	cols <- c("mistyrose",  "cornsilk", "lightblue", "darkblue", "gray", "black")

	p.dat <- p.dat[c(5,6,3,4,1,2), ];
	
	for(i in 1:dim(p.dat)[2])
	{
		bar.w <- 0.8/NROW(p.dat);
		for(j in 1:NROW(p.dat))
			rect(0.6+1*(i-1)+(j-1)*bar.w, 0, 0.6+1*(i-1)+j*bar.w, p.dat[j,i], col = cols[j], border = "black" );
	}
	axis(1, at=c(1:3), labels=c(500, 1000,2500), cex.axis=0.75)

	text( 0.4, 1.00, letter, cex=1.2, adj=c(0, 0.5), font=2);

	y.offset <- 0.025;
	#x.offset <- 2.6-bar.w*3;
	x.offset <- 2.35;
	
	if(legend)
	{
		rect( x.offset + 0.5, y.offset + 1.01, x.offset + 0.62, y.offset + 0.95, col=cols[1], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.98, expression('LSKAT'), cex=0.6, adj=c(0, 0.5));

		rect( x.offset + 0.5, y.offset + 0.93, x.offset + 0.62, y.offset + 0.87, col=cols[2], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.90, "LBT", cex=0.6, adj=c(0, 0.5));

		rect( x.offset + 0.5, y.offset + 0.85, x.offset + 0.62, y.offset + 0.79, col=cols[3], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.82,  expression('SKAT'[AVG]), cex=0.6, adj=c(0, 0.5));

		rect( x.offset + 0.5, y.offset + 0.77, x.offset + 0.62, y.offset + 0.71, col=cols[4], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.74,  expression('BT'[AVG]), cex=0.6, adj=c(0, 0.5));

		rect( x.offset + 0.5, y.offset + 0.69, x.offset + 0.62, y.offset + 0.63, col=cols[5], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.66,  expression('SKAT'[BL]), cex=0.6, adj=c(0, 0.5));

		rect( x.offset + 0.5, y.offset + 0.6, x.offset + 0.62, y.offset + 0.55, col=cols[6], border="black")
		text( x.offset + 0.7-0.02, y.offset + 0.58,  expression('BT'[BL]), cex=0.6, adj=c(0, 0.5));
	}
}

load.lskat.ret <- function(test.no, data.folder="./", loops=100)
{
	load( paste( data.folder, simu.rdata[test.no],  ".",  500, ".rdata", sep="") )
	r.power0 <- r.power6;
	load( paste( data.folder, simu.rdata[test.no],  ".", 1000, ".rdata", sep="") );
	r.power0 <- rbind( r.power0, r.power6) ;
	load( paste( data.folder, simu.rdata[test.no],  ".", 2500, ".rdata", sep="") );
	r.power0 <- rbind( r.power0, r.power6);

	pp1<- r.power0[which(r.power0[,2]==1), c(6,7,8,9,3,4)];
	pp1<- t(pp1/loops);
	pp2<- r.power0[which(r.power0[,2]==2), c(6,7,8,9,3,4)];
	pp2<- t(pp2/loops);
	pp3<- r.power0[which(r.power0[,2]==3), c(6,7,8,9,3,4)];
	pp3<- t(pp3/loops);
	
	return(list(pp1=pp1, pp2=pp2, pp3=pp3));
}

draw.3by3<-function( test.no1, test.no2, test.no3, pdf.file, data.folder, loops=100)
{
	r.power0 <- c();

	#pdf(pdf.file, width=6, height=6);
	postscript(pdf.file, width=6, height=6);
	
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
		c( 0,    10,10,10, 0 ),
		c( 0,    11,11,11, 0 ),
		c( 0,   1:3, 0 ),
		c( 0,    12,12,12, 0 ),
		c( 0,   4:6, 0 ),
		c( 0,    13,13,13, 0 ),
		c( 0,    7:9, 0 ),
		c( 0,    14,14,14, 0 ),
		c( rep(0, 5) ));
	
	nf <- layout(lay.mat, c(0.01, 2,  2,  2, 0.01), c(0.01, 0.5, 0.5, 4.5, 0.5, 4.5, 0.5, 4.5, 0.5, 0.01 ) );

	par(mar=c( 1.5, 3, 0.75, 0.5 ));

	r1 <- load.lskat.ret( test.no1, data.folder, loops=loops)
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r1$pp1, ylim=c(0,1.05),  "", ylab=T);
	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(r1$pp2, ylim=c(0,1.05),  "" );
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r1$pp3, ylim=c(0,1.05),  "", legend=T)

	par(mar=c( 1.5, 3, 0.75, 0.5 ));

	r2 <- load.lskat.ret( test.no2, data.folder, loops=loops)
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r2$pp1, ylim=c(0,1.05), "", ylab=T);
	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(r2$pp2, ylim=c(0,1.05), "" );
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r2$pp3, ylim=c(0,1.05), "" )

	par(mar=c( 1.5, 3, 0.75, 0.5 ));

	r3 <- load.lskat.ret( test.no3, data.folder, loops=loops)
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r3$pp1, ylim=c(0,1.05), "", ylab=T );
	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(r3$pp2, ylim=c(0,1.05), "" );
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r3$pp3, ylim=c(0,1.05), "" )

	x.offset <- 0.1;
	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 1.5, 0.5, "Model I", adj=c(0.5,0.5), cex=1.4, font=2);

	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 0.5, 0.25, "100% Positive", adj=c(0.5,0.5), cex=1);
	text(x.offset + 1.5, 0.25, "80% Positive",  adj=c(0.5,0.5), cex=1);
	text(x.offset + 2.5, 0.25, "50% Positive",  adj=c(0.5,0.5), cex=1);
	
	
	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F);
	text(x.offset + 1.5, 0.5, "Model II", adj=c(0.5,0.5), cex=1.4, font=2);

	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F);
	text(x.offset + 1.5, 0.5, "Model III", adj=c(0.5,0.5), cex=1.4, font=2);

	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 0.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	text(x.offset + 1.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	text(x.offset + 2.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	
	dev.off();
}


draw.2by3<-function( test.no1, test.no2, pdf.file, data.folder, loops=100)
{
	r.power0 <- c();

	#pdf(pdf.file, width=6, height=4);
	postscript(pdf.file, width=6, height=4);
	
	opar <- par(no.readonly=TRUE);
	par(pty = "m", lwd=1, oma=c(0,0,0,0) );

	lay.mat<-rbind(
		c( 0,   rep(0,3),  0 ),
		c( 0,    7,7,7, 0 ),
		c( 0,    8,8,8, 0 ),
		c( 0,   1:3, 0 ),
		c( 0,    9,9,9, 0 ),
		c( 0,   4:6, 0 ),
		c( 0,    10,10,10, 0 ),
		c( rep(0, 5) ));
	
	nf <- layout(lay.mat, c(0.01, 2,  2,  2, 0.01), c(0.01, 0.5, 0.5, 4.5, 0.5, 4.5, 0.5, 0.01 ) );

	par(mar=c( 1.5, 3, 0.75, 0.5 ));

	r1 <- load.lskat.ret( test.no1, data.folder, loops=loops)
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r1$pp1, ylim=c(0,1.05),  "", ylab=T);
	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(r1$pp2, ylim=c(0,1.05),  "" );
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r1$pp3, ylim=c(0,1.05),  "", legend=T)

	par(mar=c( 1.5, 3, 0.75, 0.5 ));

	r2 <- load.lskat.ret( test.no2, data.folder, loops=loops)
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r2$pp1, ylim=c(0,1.05), "", ylab=T);
	par(mgp=c(1.4, 0.3, 0 ), tck=-0.03);
	bar.subplot(r2$pp2, ylim=c(0,1.05), "" );
	par(mgp=c(1.4, 0.3, 0), tck=-0.03);
	bar.subplot(r2$pp3, ylim=c(0,1.05), "" )

	x.offset <- 0.1;
	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 1.5, 0.5, "Model IV", adj=c(0.5,0.5), cex=1.4, font=2);

	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 0.5, 0.25, "100% Positive", adj=c(0.5,0.5), cex=1);
	text(x.offset + 1.5, 0.25, "80% Positive",  adj=c(0.5,0.5), cex=1);
	text(x.offset + 2.5, 0.25, "50% Positive",  adj=c(0.5,0.5), cex=1);
	
	
	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F);
	text(x.offset + 1.5, 0.5, "Model V", adj=c(0.5,0.5), cex=1.4, font=2);

	par(mar=c( 0, 0, 0, 0 ));
	par(mgp=c(0,0,0));
	plot(1,1, type="n", xlim=c(0,3), ylim=c(0,1), xaxs="i", yaxs="i",axes=F );
	text(x.offset + 0.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	text(x.offset + 1.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	text(x.offset + 2.5, 0.6, "Sample Size", adj=c(0.5,0.5), cex=1);
	
	dev.off();
}

#draw.3by3( 1, 5, 4, "Power-154.pdf", "simu-100/", 100);
#draw.2by3( 14, 12, "Power-EC.pdf", "simu-100/", 100);


#draw.3by3( 1, 5, 4, "Power-154.ps", "simu-100/", 100);
#draw.2by3( 14, 12, "Power-EC.ps", "simu-100/", 100);


#draw.3by3( 1, 5, 4, "Power-154.pdf", "./", 1000);
#draw.2by3( 14, 12, "Power-EC.pdf", "./", 1000);

draw.3by3( 1, 5, 4, "Power-154.ps", "./", 1000);
draw.2by3( 14, 12, "Power-EC.ps", "./", 1000);
