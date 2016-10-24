library(LSKAT);

p0 <- longskat_gene_simulate( plink.format=F, power.test=T, n.gene=1);

r.model <- longskat_est_model( p0$phe.long, p0$phe.cov, g.maxiter=3, verbose=T);
print(r.model);

for(i in 1:NCOL(p0$snp.mat[[1]]))
{
	r.lskat <- longskat_snp_test(r.model, p0$snp.mat[[1]][,i]);
	print(r.lskat);
}

p1 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-simu-snp", power.test=T, n.gene=2);
r.lskat <- longskat_snp_plink( p1$file.plink.bed, p1$file.plink.bim, p1$file.plink.fam, p1$file.phe.long, p1$file.phe.cov, NULL, p1$file.gene.set );
print(r.lskat);


plot(r.lskat, "tmp-simu-snp.pdf", "title")