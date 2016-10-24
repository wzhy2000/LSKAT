library(LSKAT);

p0 <- longskat_gene_simulate();
r.model0 <- longskat_est_model( p0$phe.long, p0$phe.cov, g.maxiter=3, verbose=T)
print(r.model0);

r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[1]]);
print(r.lskat0);

r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[1]], rare.cutoff=0.05, test.type="Joint");
print(r.lskat0);

r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[1]], rare.cutoff=0.05, test.type="Common.Only");
print(r.lskat0);

r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[1]], rare.cutoff=0.05, test.type="Common.onlY");
print(r.lskat0);

r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[1]], rare.cutoff=0.05, test.type="rare.onlY");
print(r.lskat0);

p2 <- longskat_gene_simulate(power.test=T, geno.miss=0.05,  pheno.miss=0.08, pheno.cov="SAD", intercept=TRUE);
r.model2 <- longskat_est_model( p2$phe.long, p2$phe.cov, phe.time = NULL, g.maxiter=5, verbose=F)
print(r.model2);
if(!is.null(r.model2))
{
	r.lskat2 <- longskat_gene_test(r.model2, p2$snp.mat[[1]]);
	print(r.lskat2);

	r.lskat2 <- longskat_gene_test(r.model2, p2$snp.mat[[1]], rare.cutoff=0.05, test.type="Common.only");
	print(r.lskat2);

	r.lskat2 <- longskat_gene_test(r.model2, p2$snp.mat[[1]], rare.cutoff=0.05, test.type="Rare.only");
	print(r.lskat2);
}



p1 <- longskat_gene_simulate(power.test=F, geno.miss=0.02,  pheno.miss=0.05, pheno.cov="cs");
r.model1 <- longskat_est_model( p1$phe.long, p1$phe.cov, phe.time = NULL, g.maxiter=3, verbose=F)
print(r.model1);

if(!is.null(r.model1))
{
	r.lskat1 <- longskat_gene_test(r.model1, p1$snp.mat[[1]]);
	print(r.lskat1);


	r.lskat1 <- longskat_gene_test(r.model1, p1$snp.mat[[1]], rare.cutoff=0.05, test.type="Common.only");
	print(r.lskat1);

	r.lskat1 <- longskat_gene_test(r.model1, p1$snp.mat[[1]], rare.cutoff=0.05, test.type="Rare.only");
	print(r.lskat1);
}