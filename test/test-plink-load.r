library(LSKAT);

p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-plink-load", power.test=T, n.gene=20);

r.model <- longskat_est_model( p0$phe.long, p0$phe.cov, g.maxiter=3, verbose=T);

plk <- longskat_plink_load ( p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.gene.set,  verbose=TRUE);

snps <- longskat_get_gene( plk, 1:2, snp.impute="mean" )

str(snps);

longskat_gene_test( r.model, snps$snp.mat[[1]], verbose=T )

longskat_gene_test( r.model, snps$snp.mat[[2]], verbose=T )

