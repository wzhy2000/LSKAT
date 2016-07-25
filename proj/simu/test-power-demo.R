source("long-simu-testplan.r");

# $R  --vanilla --quite power=1 nloop=1000 phe.sample=500  phe.dist=mn phe.cov=sad phe.par=0.7,0.4 ret.rdata=power.L13.mn.sad.500.rdata sig.a=0.4 sig.e=0.4  < long-simu-testplan.r > power-L13-mn-sad-500.out
testR( 1, "power.LX.mn.sad.500.rdata",
		nloop    = 10, 
		nrep     = 1, 
		nsample  = 500, 
		phe.dist = "mn", 
		phe.cov  = "sad", 
		rho      = 0.7,
		sig_a    = 0.8, 
		sig_b    = 0.8, 
		sig_e    = 0.8,
		phe.par1 = 0.7, 
		phe.par2 = 0.8,
		ncores   = 7)
