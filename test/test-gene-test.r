library(LSKAT);

p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-test", power.test=F );
r.ml <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3, verbose=T, method="ML")
print(r.ml);

r.reml <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3,  verbose=T, method="REML")
print(r.reml);

r.reml0 <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3, intercept=T, verbose=T, method="REML", time.cov=2)

for(i in 1:length(p0$snp.mat))
{
	r.lskat0 <- longskat_gene_test(r.reml, p0$snp.mat[[i]], snp.impute="mean");
	print(r.lskat0);
}

for(i in 1:length(p0$snp.mat))
{
	r.lskat1 <- longskat_gene_test(r.reml, p0$snp.mat[[i]], snp.impute="random", weights.common=c(0.5,1), weights.rare=c(1,5), rare.cutoff=0.05, test.type="Common.Only" );
	print(r.lskat1);
}


p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-test", power.test=T );
r.ml <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3, verbose=T, method="ML")
print(r.ml);

r.reml <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3,  verbose=T, method="REML")
print(r.reml);

r.reml0 <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3, intercept=T, verbose=T, method="REML", time.cov=2)
print(r.reml0);


for(i in 1:length(p0$snp.mat))
{
	r.lskat0 <- longskat_gene_test(r.reml, p0$snp.mat[[i]], snp.impute="mean");
	print(r.lskat0);
}

for(i in 1:length(p0$snp.mat))
{
	r.lskat1 <- longskat_gene_test(r.reml, p0$snp.mat[[i]], snp.impute="random", weights.rare=c(1,5), rare.cutoff=0.10, test.type="Rare.Only" );
	print(r.lskat1);
}

