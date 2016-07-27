source("long-simu-test-rare2.r");

if(!exists("ncores")) ncores<-16;

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
    snprange    = c(10*1000, 60*1000),
	a.level     = 10^(-6),
	w.common    = c(0.5, 0.5),
	w.rare      = c(1,25),
	effect.sd   = 0.05,
	# 2015 test: rare.c1=0.08, common.c1=0.12, max.rare.causal=10, max.common.causal=4,
	rare.count.range = c(30, 32),
	common.count.range = c(14, 16),
	rare.c1     = 0.08/3,
	common.c1   = 0.12/3,
	max.rare.causal = 30,
	max.common.causal = 12,
    positive.effect = 1,
	rare.cutoff = 0.05 );

n.loop <- 800;
n.rep <- 1

#testR( 1, par, "power.L1.1k.rare.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L1.1k.rare.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L1.1k.rare.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.2
par$sig_b <- 0.8
par$sig_e <- 0.2
#testR( 1, par, "power.L4.1k.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L4.1k.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L4.1k.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.2
par$sig_e <- 0.8
#testR( 1, par, "power.L5.1k.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L5.1k.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L5.1k.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);


par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
#testR( 1, par, "power.L12.1k.mn.sad.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L12.1k.mn.sad.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L12.1k.mn.sad.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);


par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
#testR( 1, par, "power.L14.1k.mn.cm.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
#testR( 1, par, "power.L14.1k.mn.cm.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L14.1k.mn.cm.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);




if(0)
{
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
	effect.sd   = 0.05,
	# 2015 test: rare.c1=0.08, common.c1=0.12, max.rare.causal=10, max.common.causal=4,
	rare.c1     = 0.08,
	common.c1   = 0.12,
	max.rare.causal = 0,
	max.common.causal = 4,
    positive.effect = 1,
	rare.cutoff = 0.05 );

n.loop <- 200;
n.rep <- 1

testR( 1, par, "power.L1.200.common.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.common.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.common.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);


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
	effect.sd   = 0.05,
	# 2015 test: rare.c1=0.08, common.c1=0.12, max.rare.causal=10, max.common.causal=4,
	rare.c1     = 0.08,
	common.c1   = 0.12,
	max.rare.causal = 4,
	max.common.causal = 4,
    positive.effect = 1,
	rare.cutoff = 0.05 );

n.loop <- 200;
n.rep <- 1

testR( 1, par, "power.L1.200.snp4.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.snp4.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.snp4.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);


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
	effect.sd   = 0.05,
	# 2015 test: rare.c1=0.08, common.c1=0.12, max.rare.causal=10, max.common.causal=4,
	rare.c1     = 0.08,
	common.c1   = 0.12,
	max.rare.causal = 10,
	max.common.causal = 10,
    positive.effect = 1,
	rare.cutoff = 0.05 );

n.loop <- 200;
n.rep <- 1

testR( 1, par, "power.L1.200.snp10.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.snp10.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L1.200.snp10.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

}



if(0)
{
par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
testR( 1, par, "power.L2.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L2.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L2.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.2
par$sig_b <- 0.2
par$sig_e <- 0.2
testR( 1, par, "power.L3.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L3.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L3.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.2
par$sig_b <- 0.8
par$sig_e <- 0.2
testR( 1, par, "power.L4.1k.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L4.1k.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L4.1k.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.2
par$sig_e <- 0.8
testR( 1, par, "power.L5.1k.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L5.1k.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L5.1k.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
testR( 1, par, "power.L6.mt.ar1.500.rdata",  nsample =  500, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L6.mt.ar1.1000.rdata", nsample = 1000, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L6.mt.ar1.2500.rdata", nsample = 2500, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
testR( 1, par, "power.L7.mt.ar1.500.rdata",  nsample =  500, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L7.mt.ar1.1000.rdata", nsample = 1000, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L7.mt.ar1.2500.rdata", nsample = 2500, phe.dist = "mt", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
testR( 1, par, "power.L8.msn.ar1.500.rdata",  nsample =  500, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L8.msn.ar1.1000.rdata", nsample = 1000, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L8.msn.ar1.2500.rdata", nsample = 2500, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
testR( 1, par, "power.L9.msn.ar1.500.rdata",  nsample =  500, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L9.msn.ar1.1000.rdata", nsample = 1000, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L9.msn.ar1.2500.rdata", nsample = 2500, phe.dist = "msn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
par$par1 <- 0.5
par$par2 <- 0.8
par$par3 <- 0.5
par$par4 <- 0.8
par$par5 <- 0.6
testR( 1, par, "power.L10.mmn.ar1.500.rdata",  nsample =  500, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L10.mmn.ar1.1000.rdata", nsample = 1000, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L10.mmn.ar1.2500.rdata", nsample = 2500, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
par$par1 <- 0.5
par$par2 <- 0.8
par$par3 <- 0.3
par$par4 <- 0.4
par$par5 <- 0.7
testR( 1, par, "power.L11.mmn.ar1.500.rdata",  nsample =  500, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L11.mmn.ar1.1000.rdata", nsample = 1000, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L11.mmn.ar1.2500.rdata", nsample = 2500, phe.dist = "mmn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
testR( 1, par, "power.L12.1k.mn.sad.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L12.1k.mn.sad.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L12.1k.mn.sad.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
testR( 1, par, "power.L13.mn.sad.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L13.mn.sad.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L13.mn.sad.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "sad", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.8
par$sig_b <- 0.8
par$sig_e <- 0.8
testR( 1, par, "power.L14.1k.mn.cm.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L14.1k.mn.cm.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L14.1k.mn.cm.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.4
par$sig_b <- 0.4
par$sig_e <- 0.4
testR( 1, par, "power.L15.mn.cm.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L15.mn.cm.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L15.mn.cm.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "cm", nloop = n.loop, nrep = n.rep, ncores = ncores);

par$sig_a <- 0.1
par$sig_b <- 1.48
par$sig_e <- 1.67
par$rho   <- 0.87

testR( 1, par, "power.L16.mn.ar1.500.rdata",  nsample =  500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L16.mn.ar1.1000.rdata", nsample = 1000, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);
testR( 1, par, "power.L16.mn.ar1.2500.rdata", nsample = 2500, phe.dist = "mn", phe.cov  = "ar1", nloop = n.loop, nrep = n.rep, ncores = ncores);

}
