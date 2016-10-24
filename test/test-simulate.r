library(LSKAT);

p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-plink-simulate", power.test=T );
r.lskat1 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, 
								options=list(g.maxiter=3, plink.path="plink"));

r.lskat1


p0 <- longskat_gene_simulate( plink.format=T, file.plink.prefix="tmp-plink-simulate", power.test=F );
r.lskat2 <- longskat_gene_plink(p0$file.plink.bed, p0$file.plink.bim, p0$file.plink.fam, p0$file.phe.long, p0$file.phe.cov, NULL, p0$file.gene.set, 
								options=list(g.maxiter=3, plink.path="plink"));


p1 <- longskat_gene_simulate( plink.format=F, power.test=F );
p2 <- longskat_gene_simulate( plink.format=F, power.test=F,
	n.minsect = 3000,
	n.maxsect = 3000,
	n.sample = 1000,
	n.time = 8,
	n.gene = 20,
	geno.miss = 0.00,
	pheno.miss = 0.00,
	pheno.dist = "MT",
	pheno.cov = "SAD",
	intercept = T);

p3 <- longskat_gene_simulate( plink.format=F, power.test=T,
	n.minsect = 3000,
	n.maxsect = 3000,
	n.sample = 1000,
	n.time = 8,
	n.gene = 1,
	geno.miss = 0.01,
	pheno.miss = 0.01,
	pheno.dist = "MT",
	pheno.cov = "CS",
	intercept = T);

p4 <- longskat_gene_simulate( plink.format=F, power.test=T,
	n.minsect = 3000,
	n.maxsect = 3000,
	n.sample = 1000,
	n.time = 8,
	n.gene = 1,
	geno.miss = 0.01,
	pheno.miss = 0.01,
	pheno.dist = "MSN",
	pheno.cov = "AR1",
	intercept = T);
	

p5 <- longskat_gene_simulate( plink.format=F, power.test=T,
	n.minsect = 3000,
	n.maxsect = 3000,
	n.sample = 1000,
	n.time = 8,
	n.gene = 1,
	geno.miss = 0.01,
	pheno.miss = 0.01,
	pheno.dist = "MMN",
	pheno.cov = "AR1",
	intercept = T);
	
p6 <- longskat_gene_simulate( plink.format=F, power.test=T,
	n.minsect = 3000,
	n.maxsect = 3000,
	n.sample = 1000,
	n.time = 8,
	n.gene = 1,
	geno.miss = 0.01,
	pheno.miss = 0.01,
	pheno.dist = "MSN",
	pheno.cov = "SAD",
	intercept = T);

p7 <- longskat_gene_simulate( plink.format=F, power.test=T,
	n.minsect = 10000,
	n.maxsect = 30000,
	n.sample = 1000,
	n.time = 12,
	n.gene = 10,
	geno.miss = 0.00,
	pheno.miss = 0.00,
	pheno.dist = "MN",
	pheno.cov = "AR1",
	intercept = T);
	