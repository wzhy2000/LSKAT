read_phe_cov<-function( file.phe.long, file.phe.time, file.phe.cov, PF.gen)
{
	phe.long <- read.csv(file.phe.long, header=T, stringsAsFactors=F, row.names=1);
	cat("  PHE LONG =", file.phe.long, "\n");
	cat("* Individuals =", NROW(phe.long), "\n");
	cat("* Times =", NCOL(phe.long), "\n");
	cat("* Mean =",  colMeans(phe.long, na.rm=T), "\n");
	cat("* SD =",    colSds(phe.long, na.rm=T),"\n");
	idx.na <- which( rowSums(is.na(phe.long)) == NCOL(phe.long) );
	if( length(idx.na)>0) phe.long <- phe.long[ -idx.na, ];

	phe.time <- NULL;
	if (!is.null(file.phe.time))
	{
		phe.time <- read.csv(file.phe.time, header=T, stringsAsFactors=F, row.names=1);
		cat("  PHE TIME =", file.phe.time, "\n");
		cat("* Individuals =", NROW(phe.time), "\n");
		cat("* Times =", NCOL(phe.time), "\n");
		cat("* Mean =",  colMeans(phe.time, na.rm=T), "\n");
		cat("* SD =",    colSds(phe.time, na.rm=T),"\n");
		idx.na <- which( rowSums( is.na(phe.time))==NCOL(phe.time) );
		if( length(idx.na)>0) phe.time <- phe.time[ -idx.na, ];
	}

	phe.cov <- read.csv(file.phe.cov, header=T, stringsAsFactors=F, row.names=1);
	cat("  PHE COV =", file.phe.cov, "\n");
	cat("* Individuals =", NROW(phe.cov), "\n");
	cat("* Covariate =", NCOL(phe.cov), "\n");
	cat("* Mean =",  colMeans(phe.cov, na.rm=T), "\n");
	cat("* SD =",    colSds(phe.cov, na.rm=T), "\n");

	ids.phe <- intersect(rownames(phe.long), rownames(phe.cov) );
	if(!is.null(phe.time))
		ids.phe <- intersect(ids.phe, rownames(phe.time) );

	ids.set <- intersect(ids.phe, get_gen_individuals(PF.gen) );
	cat("  COMMON Individuals=", length(ids.set), "\n");

	#eg. c(10:1)[match(c(4, 6,8,2,3), c(10:1))]

	idx.long <- match( ids.set, rownames(phe.long) );
	phe.long <- phe.long[idx.long, ,drop=F];

	idx.cov <- match( ids.set, rownames(phe.cov) );
	phe.cov <- phe.cov[idx.cov, , drop=F];

	if(!is.null(phe.time))
	{
		idx.time <- match( ids.set, rownames(phe.time) );
		phe.time <- phe.time[idx.time, ,drop=F];
	}

	sync_gen_individuals(PF.gen, ids.set)

	if( !is.null(phe.time) && !all( rownames(phe.long) == rownames(phe.time) ) )
		stop("! ID MATCH ERROR between PHE.LONG and PHE.TIME. \n");

	if (!( all( rownames(phe.long)==rownames(phe.cov)) && all( rownames(phe.long)== get_gen_individuals(PF.gen) ) ) )
		stop("! ID MATCH ERROR among 3 files( PHE.LONG, PHE.COV, PLINK.FAM). \n");

	return(list( phe.long=phe.long, phe.time=phe.time, phe.cov = phe.cov));
}