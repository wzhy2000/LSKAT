library(LSKAT);

p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-plink", power.test=F, n.gene=20);
r.lskat0 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, gene.set=c(5:10), verbose=T );
r.lskat0;


p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-gene-plink", power.test=T );
r.model0 <- longskat_est_model( p0$phe.long, p0$phe.cov, phe.time = NULL, g.maxiter=3, verbose=T, method="ML")
print(r.model0);
for(i in 1:10)
{
	r.lskat0 <- longskat_gene_test(r.model0, p0$snp.mat[[i]], snp.impute="mean");
	print(r.lskat0);

}

r.lskat1 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, 
								options=list(g.maxiter=3, plink.path="plink", est.method = "ML"), verbose=T );
r.lskat1;


df <- summary(r.lskat1);
plot( r.lskat1);
plot( r.lskat1, pdf.file="tmp-plink-lskat1.pdf", title="TEST PLINK 1", bonferroni=F);

r.lskat2 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, gene.set=c(1:5),
							    options=list(time.cov=2, g.maxiter=6, rare.cutoff=0.05, snp.impute="random", intercept=TRUE));
r.lskat2;
df <- summary(r.lskat2);
plot( r.lskat2);
plot( r.lskat2, pdf.file="tmp-plink-lskat2.pdf", title="TEST PLINK 2", bonferroni=T);

r.lskat3<-longskat_snp_plink( p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, file.gene.set=p0$file.gene.set, 
				    snp.set=c(1:199),
				    options=list(time.cov=1, g.maxiter=6, n.cpu=2, test.type="Rare.only" ));
r.lskat3;
df <- summary(r.lskat3);
plot( r.lskat3);
plot( r.lskat3, pdf.file="tmp-plink-lskat3.pdf", title="TEST PLINK 3", bonferroni=T);

unlink( c(p$file.plink.bed, p$file.plink.bim, p$file.plink.fam, p$file.gene.set, p$file.phe.long, p$file.phe.cov ) ); 