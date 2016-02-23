library(LSKAT);

g.snp.hap1 <- "/ycga-ba/home/zw224/f/R/LSKAT/proj/simu/skat-test-1.hap"
g.snp.pos1 <- "/ycga-ba/home/zw224/f/R/LSKAT/proj/simu/skat-test-1.pos"

p <- longskat_gene_simulate(g.snp.hap1, g.snp.pos1, plink.format=T, file.plink.prefix="test-plink-lskat");

r.lskat1 <- longskat_gene_plink(p$file.plink.bed, p$file.plink.bim, p$file.plink.fam, p$file.gene.set, p$file.phe.long, p$file.phe.cov, 
								options=list(g.maxiter=6));
r.lskat1;
summary(r.lskat1);
plot( r.lskat1);
plot( r.lskat1, pdf.file="test-plink-lskat1.pdf", title="TEST PLINK 1", bonferroni=F);

r.lskat2 <- longskat_gene_plink(p$file.plink.bed, p$file.plink.bim, p$file.plink.fam, p$file.gene.set, p$file.phe.long, p$file.phe.cov, gene.range=c(1:5)
							    options=list(y.cov.time=2, g.maxiter=6));
r.lskat2;
summary(r.lskat2);
plot( r.lskat2);
plot( r.lskat2, pdf.file="test-plink-lskat2.pdf", title="TEST PLINK 2", bonferroni=T);

r.lskat3<-longskat_snp_plink( p$file.plink.bed, p$file.plink.bim, p$file.plink.fam, p$file.gene.set, p$file.phe.long, p$file.phe.cov, 
				    snp.range=c(1:199), 
				    options=list(y.cov.time=2, g.maxiter=6));
r.lskat3;
summary(r.lskat3);
plot( r.lskat3);
plot( r.lskat3, pdf.file="test-plink-lskat3.pdf", title="TEST PLINK 3", bonferroni=T);

unlink( c(p$file.plink.bed, p$file.plink.bim, p$file.plink.fam, p$file.gene.set, p$file.phe.long, p$file.phe.cov ) );