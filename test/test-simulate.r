library(LSKAT);

g.snp.hap1 <- "/ycga-ba/home/zw224/f/R/LSKAT/proj/simu/skat-test-1.hap"
g.snp.pos1 <- "/ycga-ba/home/zw224/f/R/LSKAT/proj/simu/skat-test-1.pos"

p <- longskat_gene_simulate(g.snp.hap1, g.snp.pos1);

r.model <- longskat_est_model( p$phe.y, p$phe.cov, y.time = NULL, g.maxiter=3, debug=T)
print(r.model);

r.lskat <- longskat_gene_run(r.model, p$snp);
print(r.lskat);

save(p, r.model, r.lskat, file="test-simulate.rdata");

library(LSKAT);
load("test-simulate.rdata");

for(i in 1:NCOL(p$snp))
{
	r.lskat <- longskat_snp_run(r.model, p$snp[,i]);
	print(r.lskat);
}
