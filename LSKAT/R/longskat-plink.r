setRefClass("PLINK.refer",
	fields = list(
			options        = "list",
			snp            = "list",
			gen.list       = "list",
			ind.list       = "list"),
			
	methods = list(
     	
		show = function() 
		{
       		cat("Reference PLINK class", classLabel(class(.self)), "\n")

			cat("PLINK BED=", options$file.plink.bed,"\n")
			cat("PLINK BIM=", options$file.plink.bim,"\n")
			cat("PLINK FAM=", options$file.plink.fam,"\n")
			cat("PLINK Path=", options$plink.path,"\n")
			cat("Individual=", NROW( snp$fam ),"\n")
			cat("Gene count=", gen.list$len,"\n")
			cat("SNP count=",  NROW( snp$bim),"\n")

	    }
	)    
);

snp_impute<-function(snp.mat, impute="mean")
{
	snp.imp <- snp.mat;

	for(i in 1:NCOL(snp.mat) )
	{
		s.mat.i <- snp.mat[,i] ;
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			if(impute=="mean")
			{
				s.mat.i[s.miss] <- mean(s.mat.i, na.rm=T);
			}
			else
			{
				n.s0 <- length( which( s.mat.i == 0 ) );
				n.s1 <- length( which( s.mat.i == 1 ) );
				n.s2 <- length( which( s.mat.i == 2 ) );
				n.s  <- length(s.mat.i)

				r.miss<- runif( length(s.miss) );
				r.snp <- rep(2, length(s.miss));
				r.snp[r.miss <= n.s0/n.s ]<-0;
				r.snp[r.miss <= (n.s0 + n.s1)/n.s ]<-1;
				s.mat.i[s.miss] <- r.snp;
			}
		}

		if (mean(s.mat.i)/2>0.5) s.mat.i <- 2 - s.mat.i;

		snp.imp[,i] <- s.mat.i;
	}

	return(snp.imp);
}

colSds<-function(mat, na.rm=T)
{
	r<-c();
	for(i in 1:dim(mat)[2])
		r <- c(r, sd(mat[,i], na.rm=na.rm));
	return(r);
}

load_gene_plink <- function(file.plink.bed, file.plink.bim, file.plink.fam, individuals, snps, plink)
{
	if( !is.null(plink) && !is.na(plink) )
	{
		tmp <- tempfile(pattern = "LSKAT.temp.");
		tmp.file.ind <- paste(tmp, "ind", sep=".");
		tmp.file.snp <- paste(tmp, "snp", sep=".");
		tmp.all <- paste(tmp, "*", sep=".");

		tb.fam <- read.table(file.plink.fam, header=F);
		idx.indi <- match(individuals, tb.fam[,2]);
		write.table( tb.fam[idx.indi, c(1,2)], file=tmp.file.ind, col.names=F, row.names=F, quote=F);
		write.table( snps, file=tmp.file.snp, col.names=F, row.names=F, quote=F);

		plink.cmd <- paste(plink, "--bed", file.plink.bed, "--fam", file.plink.fam, "--bim", file.plink.bim, "--extract", tmp.file.snp, "--keep", tmp.file.ind, "--out", tmp, "--make-bed");
		plink.status <- system(plink.cmd, wait=T, intern=T);

		snp.mat <- read.plink( paste( tmp, "bed", sep="."), paste(tmp, "bim", sep="."), paste(tmp, "fam", sep=".") );
		unlink(c(tmp.all, tmp.file.ind, tmp.file.snp));
	}
	else
	{
		snp.mat <- read.plink( file.plink.bed,  file.plink.bim, file.plink.fam);
		idx.fam <- match( individuals, snp.mat$fam$member );
		snp.mat$genotypes<- snp.mat$genotypes[idx.fam,]
		snp.mat$fam      <- snp.mat$fam[idx.fam,]
	}

	return(snp.mat);
}


get_gen_group<-function(gen.list, idx)
{
	gen.name <- gen.list$genes[idx];
	snp.idx <- which(gen.list$snps[,1]==gen.name);
	return(list(name=gen.name, snps=gen.list$snps[snp.idx,2]))
}

get_gen_family<-function(gen.lib, gen.name)
{
	snp.idx <- which(gen.lib$snps[,1]==gen.name);
	if (length(snp.idx)==0)
		return(NULL)
	else
		return(list(name=gen.name, snps=gen.lib$snps[snp.idx,2]));
}

get_gen_individuals<-function(PF.gen)
{
	return( as.character(PF.gen$ind.list$member[,2]) );
}

sync_gen_individuals<-function(PF.gen, ids.set)
{
	cat("* PLINK (", NROW(PF.gen$ind.list$member) - length(ids.set), ") individuals are removed.\n");

	idx.fam <- match( ids.set, PF.gen$ind.list$member[,2] );
	if(!is.null(PF.gen$snp$matrix))
	{
		PF.gen$snp$matrix$genotypes<- PF.gen$snp$matrix$genotypes[idx.fam,]
		PF.gen$snp$matrix$fam      <- PF.gen$snp$matrix$fam[idx.fam,]
		PF.gen$ind.list$removed    <- setdiff(PF.gen$ind.list$member[,2], PF.gen$snp$matrix$fam );
		PF.gen$ind.list$member     <- PF.gen$ind.list$member[idx.fam, ];
	}
	else
	{
		PF.gen$ind.list$removed    <- setdiff(PF.gen$ind.list$member[,2], ids.set );
		PF.gen$ind.list$member     <- PF.gen$ind.list$member[idx.fam, ];
	}
}

get_gen_mat<-function( PF.gen, idx, impute="mean" )
{
	get_plink_mat<-function(plink, snps, impute)
	{
		snp.idx <- match(as.character(snps), as.character(plink$map[,2]));
		if (length(which(is.na(snp.idx)))>0)
			snp.idx <- snp.idx[-which(is.na(snp.idx))];

		if(length(snp.idx)==0) return(NULL);

		map <- plink$map[snp.idx, ,drop=F];
		plink.org <- as( plink$genotypes[, snp.idx, drop=F ], "numeric");

		nmiss <- apply(plink.org, 2, function(snp){sum(is.na(snp))});
		snp.imp <- snp_impute(plink.org , impute=impute)

		return(list(maf=colMeans(snp.imp)/2, snp=snp.imp, nmiss=nmiss, info=map[,c(2,1,4)]) );
	}

	gen.name <- PF.gen$gen.list$names[idx];
	snps_finding <- unique(PF.gen$gen.list$snps[which(PF.gen$gen.list$snps[,1] == gen.name), 2] );
	snp.mat <- NULL;
	
	if(!is.null(PF.gen$snp$matrix))
	{
		snp.mat <- get_plink_mat(PF.gen$snp$matrix, snps_finding, impute)
	}

	if(is.null(snp.mat))
	{
		idx.range <- c(idx-50, idx+50);
		if (idx.range[1]<1) idx.range[1] <- 1
		if (idx.range[2]>PF.gen$gen.list$len) idx.range[2] <- PF.gen$gen.list$len;

		gen.names <- PF.gen$gen.list$names[idx.range[1]:idx.range[2]];
		snps <- PF.gen$gen.list$snps[which(PF.gen$gen.list$snps[,1] %in% gen.names), 2]

		PF.gen$snp$matrix <- load_gene_plink( PF.gen$options$file.plink.bed,
				PF.gen$options$file.plink.bim,
				PF.gen$options$file.plink.fam,
				PF.gen$ind.list$member[,2],
				unique(snps),
				PF.gen$options$plink );

		snp.mat <- get_plink_mat(PF.gen$snp$matrix, snps_finding, impute)
	}

	if(!is.null(snp.mat))
		snp.mat$name <- gen.name;

	return(snp.mat);
}

get_snp_mat <- function(PF.gen, idx, impute="mean" )
{
	get_plink_snp<-function(plink, snp.name, impute)
	{
		snp.idx <- match(as.character(snp.name), as.character(plink$map[,2]));
		if (length(which(is.na(snp.idx)))>0)
			snp.name <- snp.name[-which(is.na(snp.idx))];
		
		if(length(snp.name)==0) return(NULL);

		plink.org <- as( plink$genotypes[, snp.idx, drop=F ], "numeric");
		nmiss <- apply(plink.org, 2, function(snp){sum(is.na(snp))});
		snp.imp <- snp_impute(plink.org , impute=impute)
		map <- plink$map[snp.idx, ,drop=F];

		gene.name <- "";
		if (!is.null(PF.gen$gen.list$snps))
		{
			gen.idx <- match( snp.name, PF.gen$gen.list$snps[,2])
			if (length(gen.idx)>0)
				gene.name <-  PF.gen$gen.list$snps[gen.idx[1],1];
		}
		return(list(snp=snp.imp, 
			name=snp.name, chr=map[1], loc=map[4], gene=gene.name, 
			maf=colMeans(snp.imp)/2, nmiss=nmiss, info=map[,c(2,1,4)]) );
	}
	
	snp.name <- PF.gen$snp$bim[idx,2];
	snp.mat <- NULL;
	
	if(!is.null(PF.gen$snp$matrix)  )
	{
		if( !is.na(match(snp.name, PF.gen$snp$matrix$map$snp.name ) ) )
			snp.mat <- get_plink_snp(PF.gen$snp$matrix, snp.name, impute)
	}

	if(is.null(snp.mat))
	{
		idx.range <- c(idx - 5000, idx + 5000);
		if (idx.range[1] < 1 ) idx.range[1] <- 1
		if (idx.range[2] > NROW(PF.gen$snp$bim)) idx.range[2] <- NROW(PF.gen$snp$bim);

		snp.names <- PF.gen$snp$bim[idx.range, 2];

		PF.gen$snp$matrix <- load_gene_plink( PF.gen$options$file.plink.bed,
				PF.gen$options$file.plink.bim,
				PF.gen$options$file.plink.fam,
				PF.gen$ind.list$member[,2],
				snp.names,
				PF.gen$options$plink );

		snp.mat <- get_plink_snp(PF.gen$snp$matrix, snp.name, impute)
	}

	return(	snp.mat)
}


read_gen_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, plink.path)
{
	gen.list <- list();
	if(!is.null(file.gene.set))
	{
		tb.gen <- read.table(file.gene.set, header=F, stringsAsFactors=F);
		gen.names <- unique(tb.gen[,1]);
		gen.list <- list( len=NROW(gen.names), names=gen.names, snps=tb.gen);
	}
	
	tb.fam <- read.table(file.plink.fam, header=F, stringsAsFactors=F);
	ind.list <- list( member=tb.fam[,c(1,2)], removed=c() )

	tb.bim <- read.table(file.plink.bim, header=F, stringsAsFactors=F);
 	
 	#n.snp <- get_large_file_lines( file.plink.bim);
	snp <- list()
	snp$fam <- as.data.frame(tb.fam);
	snp$bim <- as.data.frame(tb.bim);
	
	n.snp <- NROW(snp$bim);
	options <- list( plink.path=plink.path, file.plink.bed=file.plink.bed,  file.plink.bim=file.plink.bim, file.plink.fam=file.plink.fam );

	if( n.snp * 1.0 * NROW(tb.fam) < 50*1000*2000 )
	{
		snp$matrix <- snpStats::read.plink( file.plink.bed,  file.plink.bim, file.plink.fam );
	}
	else
	{
		snp$matrix <- NULL;
	}
	PLINK.refer <- getRefClass("PLINK.refer");
	PF.gen <- PLINK.refer(gen.list=gen.list, ind.list=ind.list, snp=snp, options=options);

	return(PF.gen);
}


clone_plink_refer<-function(PF.gen)
{
	PLINK.refer <- getRefClass("PLINK.refer");
	PF.gen.clone <- PLINK.refer(
				gen.list=PF.gen$gen.list,
				ind.list=PF.gen$ind.list,
				snp=list(fam=PF.gen$snp$fam, bim=PF.gen$snp$bim),
				options=PF.gen$options);
	return(PF.gen.clone);
}


#TO REMOVE
read_gen_dataset<-function( file.set, file.bim )
{
	# V2: snp
	tb.bim <- read.table(file.bim);

	# V2: snp
	tb.gen <- read.table(file.set, sep=" ", header=F);

	idx.tb <- match( as.character(tb.gen$V2), as.character(tb.bim$V2) )
	idx.gen <- c(1:NROW(tb.gen)) [ !is.na(idx.tb) ]

	genes <- unique(tb.gen[idx.gen,1]);

	return(list(len=length(genes), genes=genes, snps=tb.gen[idx.gen,]));
}


#TO REMOVE
read_gen_phe_cov<-function(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov)
{
	phe.long <- read.csv(file.phe.long, header=T, stringsAsFactors=F, row.names=1);
	idx.na <- which( rowSums(is.na(phe.long)) == NCOL(phe.long) );
	if( length(idx.na)>0) phe.long <- phe.long[ -idx.na, ];

	phe.time <- NULL;
	if (!is.null(file.phe.time))
	{
		phe.time <- read.csv(file.phe.time, header=T, stringsAsFactors=F, row.names=1);
		idx.na <- which( rowSums( is.na(phe.time))==NCOL(phe.time) );
		if( length(idx.na)>0) phe.time <- phe.time[ -idx.na, ];
	}

	phe.cov <- read.csv(file.phe.cov, header=T, stringsAsFactors=F, row.names=1);

	tb.fam <- read.table(file.plink.fam, header=F);
	ids.fam <- as.character(tb.fam[,2]);

	ids.phe <- intersect(rownames(phe.long), rownames(phe.cov) );
	if(!is.null(phe.time))
		ids.phe <- intersect(ids.phe, rownames(phe.time) );

	ids.set <- intersect(ids.phe, ids.fam);
	cat("  COMMON Individuals=", length(ids.set), "\n");

	#eg. c(10:1)[match(c(4, 6,8,2,3), c(10:1))]

	idx.long <- match( ids.set, rownames(phe.long) );
	phe.long <- phe.long[idx.long, ];

	idx.cov <- match( ids.set, rownames(phe.cov) );
	phe.cov <- phe.cov[idx.cov, ];

	if(!is.null(phe.time))
	{
		idx.time <- match( ids.set, rownames(phe.time) );
		phe.time <- phe.time[idx.time, ];
	}

	if (!all(ids.set==ids.fam) )
	{
		idx.fam <- idx.fam[ match( ids.set, ids.fam ) ];

		cat("* PLINK (", length(ids.fam) - length(ids.set), ") individuals are removed.\n");
	}

	if( !is.null(phe.time) && !all( rownames(phe.long) == rownames(phe.time) ) )
		stop("! ID MATCH ERROR between PHE.LONG and PHE.TIME. \n");

	if (!( all( rownames(phe.long)==rownames(phe.cov)) && all( rownames(phe.long)==ids.fam) ) )
		stop("! ID MATCH ERROR among 3 files( PHE.LONG, PHE.COV, PLINK.FAM). \n");

	return(list( phe.long=phe.long, phe.time=phe.time, phe.cov = phe.cov, member=idx.fam));
}

#TO REMOVE
shrink_snpmat<-function(snp.mat, gen.list, gene.range )
{
	snp.mat0 <- snp.mat;

	snp.idx <- which(!is.na(match(gen.list$snps[,1], gen.list$genes[gene.range])))
	snp.name <- unique( gen.list$snps[snp.idx,2] );

	snp.idx0 <- match( as.character(snp.name), as.character(snp.mat$map[,2]));
	if (length(which(is.na(snp.idx0)))>0)
		snp.idx0 <- snp.idx0[-which(is.na(snp.idx0))];

	if(length(snp.idx0)==0) return(NULL);

	snp.mat0$genotypes <- snp.mat$genotypes[, snp.idx0, drop=F ];
	snp.mat0$map <- snp.mat$map[snp.idx0,];

	return( snp.mat0 );
}

#public
longskat_plink_load <- function( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, plink.path=NULL, verbose=FALSE)
{
	chk.plink <- check_plink_file( file.plink.bed, file.plink.bim, file.plink.fam )
	if ( !chk.plink$bSuccess )
		stop("PLINK file can not be loaded by the snpStats package.")

	cat( "Starting to load all data files......\n");

	PF.gen <- read_gen_plink ( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, plink.path );

	return(PF.gen);
}

#public
longskat_get_gene <- function( gen.obj, gene.set, snp.impute="mean", verbose = FALSE )
{
	gene.name <- list();
	snp.list <- list();
	nmiss <- list();
	maf <- list();

	for(i in 1:length(gene.set))
	{
		gen <- try( get_gen_mat( gen.obj, gene.set[i], snp.impute ) );
		if( is.null(gen) || class(gen)=="try-error" || length(gen$maf)==0 )
		{
			if (verbose) cat("! No SNPS for Gene[", i, "]=", i, "\n");
			snp.list[[i]] <- NA;
			maf[[i]] <- NA;
			nmiss[[i]] <- NA;
			gene.name[[i]] <- NA;
		}
		else
		{
			if (verbose) cat("  Finding", NCOL(gen$snp), "SNPs...\n");
			snp.list[[i]] <- gen$snp;
			maf[[i]] <- gen$maf;
			nmiss[[i]] <- gen$nmiss;
			gene.name[[i]] <- gen$name;
		}
	}

	return(list(snp.mat=snp.list, maf=maf, nmiss=nmiss, gene.name=gene.name ));
}

