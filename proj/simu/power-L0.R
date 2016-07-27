source("long-simu-test-rare2.r");

ncores<-16;

par <- list(
	test.type   = "Joint",
    mu          = 1,
	a           = 0.5, 
	b           = 0.5, 
	rho         = 0.75,
	sig_a       = 0.8, 
	sig_b       = 0.8, 
	sig_e       = 0.8, 
	par1        = 0.7,
	par2        = 0.8,
	times       = 8, 
	intercept   = F,
	y.cov.time  = 0,
	par_t       = c(0.2, -0.08),
    snprange    = c(5*1000, 30*1000),
	a.level     = 10^(-6), 
	w.common    = c(0.5, 0.5),
	w.rare      = c(1,25),
	common.c1   = 0.12,
	rare.c1     = 0.08, 
	effect.sd   = 0.05,
	common.cause = 3,
    rare.effect = 1,
	rare.cutoff = 0.05 );	

par$sig_b <- 0;
par$rho <- 0;
n.rep <- 1
n.loop <- 200;
 
testR( 1, par, "power.L0.1k.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L0.1k.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L0.1k.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

