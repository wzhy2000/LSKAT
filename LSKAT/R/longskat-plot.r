draw_manhattan<-function( res, map.title="", sig.thres, dot.cex, y.max=NA)
{
	#par( xaxs="r",  yaxs="r");  # Extend axis limits by 4%
	#par( tck = -0.008 );
	#par( mgp = c(2,0.2,0) );
	#par( mar=c( 4, 4, 1.5, 1.5)+0.1);

	pvalues <- -log10(res[,2]);
	nrow    <- dim(res)[1];
	log10.max <- round(max(pvalues, na.rm=T))+1;

	if (is.na(y.max))
		y.max <- log10.max;

	if (length(which(pvalues>y.max))>0)
		pvalues[which(pvalues>y.max)] <- y.max;

	plot( 1,1, type="n", xlab="Chromosome", ylab=expression(-log[10](italic(p))),
		  cex.axis=0.7, xlim=c(1, nrow), ylim=c(0,  y.max), xaxt="n", main=map.title );

	p.lab <-  - c( log10( sig.thres) );

	abline( h= c(p.lab), col="gray", lwd=1, lty="dashed");

	#text( x=0, p.lab[1] + 0.1, "p=0.01", cex=0.8, srt=90, adj=c(0.5, -1)); 

	cols <- c( "darkgreen","black",  "orange",  "red", "blue", "purple");
	
	p.cex <- rep(0.5*dot.cex, length(pvalues));
	p.cex[which(pvalues>0.4*y.max)] <- 0.5*dot.cex + (pvalues[which(pvalues>y.max*0.4)]-0.4*y.max)/y.max*dot.cex;
	points( pvalues, pch=20, col=cols[ (res[,1]%%6+1)], cex=p.cex);

	x.off <- 0;
	x <- c();
	x.ps <- c();
	for(i in 1:18)
	{
		x.ps<-c( x.ps, length(which(res[,1]==i) ) );
		x <- c(x, x.ps[i]/2 + x.off);
		x.off <- x.off + x.ps[i];
	}

	x.ps<-c( x.ps, length(which(res[,1]>18 ) )  )
	x <- c(x, x.ps[19]*2/3 + x.off);
	
	axis(side=1, at=x, col="black", labels=c(paste("", 1:18, sep=""),"..."), col.axis="black", col.lab="black", cex.axis=0.4, padj=-0.2 )
}	

