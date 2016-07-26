setRefClass("PLINK.refer",
	fields = list(
			options        = "list",
			snp            = "list",
			gen.list       = "list",
			ind.list       = "list")
);

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

get_snp_plink_info<-function(idx, snp.mat, gen.tb=NA)
{
	s <- snp.mat$genotypes[, idx, drop=F ]

	s.mat <- as(s, "numeric");
	s.mat.i <- s.mat[,1];
	s.miss <- which( is.na(s.mat.i) );
	s0 <- which( s.mat.i == 0 );
	s1 <- which( s.mat.i == 1 );
	s2 <- which( s.mat.i == 2 );

	if (length(s.miss)>0)
		s.mat.i <- s.mat.i[-s.miss];

	if ( mean(s.mat.i) > 1 ) s.mat.i <- 2 -s.mat.i;
	snp.imp <- as.matrix( s.mat.i, dim=c(length(s.mat.i),1) );
	snp.maf <- sum(snp.imp)/(length(s.mat.i)*2);

	snp.map <- snp.mat$map[idx,]
	gene.name <- "";
	if (!all(is.na(gen.tb)))
	{
		gen.idx <- which( gen.tb[,2]==snp.map$snp.name)
		if (length(gen.idx)>0)
			gene.name <-  gen.tb[gen.idx[1],1];
	}

	return(list(snp=snp.imp, maf=snp.maf, name=snp.map$snp.name, chr=snp.map$chromosome, loc=snp.map$position, gene=gene.name, nmiss=length(s.miss), miss=s.miss ) );
}

get_snpstat_mat<-function(snp.mat, snps, snp.impute="mean")
{
	snps <- match(as.character(snps), as.character(snp.mat$map[,2]));
	if (length(which(is.na(snps)))>0)
		snps <- snps[-which(is.na(snps))];

	if(length(snps)==0) return(NULL);

	s.mat <- as( snp.mat$genotypes[, snps, drop=F ], "numeric");
	snp.imp <-c();
	snp.maf <- c();
	snp.names <- c();

	for(i in 1:dim(s.mat)[2])
	{
		s.mat.i <- s.mat[,i] ;
		s.miss <- which( is.na(s.mat.i) );

		if (length(s.miss)>0)
		{
			if(snp.impute=="mean")
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

		snp.imp <- rbind( snp.imp, s.mat.i );
		snp.maf <- c(snp.maf, mean(s.mat.i)/2);
		snp.names <- c(snp.names, snps[i]);
	}

	rownames(snp.imp) <- snp.names;

	map <- snp.mat$map[snps, ,drop=F];

	return(list(maf=snp.maf, snp=snp.imp, info=map[,c(2,1,4)]) );
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
	if( !is.null(plink) )
	{
		tmp.file.snp <- tempfile();
		tmp.file.ind <- tempfile();
		tmp.file.bed <- tempfile();

		tb.fam <- read.table(file.plink.fam, header=F);
		idx.indi <- match(individuals, tb.fam[,2]);
		write.table( tb.fam[idx.indi, c(1,2)], file=tmp.file.ind, col.names=F, row.names=F, quote=F);

		write.table( snps, file=tmp.file.snp, col.names=F, row.names=F, quote=F);

		cmd.plink <- paste(plink, "--bed", file.plink.bed, "--fam", file.plink.fam, "--bim", file.plink.bim, "--extract", tmp.file.snp, "--keep", tmp.file.ind, "--out", tmp.file.bed, "-make-bed");
		system( cmd.plink, internal=T, wait=T );

		snp.mat <- read.plink( paste(tmp.file.bed, "bed", sep="."), paste(tmp.file.bed, "bim", sep="."), paste(tmp.file.bed, "fam", sep=".") );
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


read_gen_phe_cov<-function(file.plink.bed, file.plink.bim, file.plink.fam, file.phe.long, file.phe.time, file.phe.cov)
{
	library(snpStats);

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

	return(list(snp.mat=snp.mat, phe.long=phe.long, phe.time=phe.time, phe.cov = phe.cov, member=idx.fam));
}

convert_simpe_to_plink <- function( snp.mat, snp.file.base )
{
	chromosome <- snp.mat[,1];
	position <- snp.mat[,2];

	# PLINK raw data: 1/2/3==> AA,AB,BB, 0==>NA
	snp.mat <- snp.mat[,-c(1,2),drop=F] + 1;

	sub.name <- colnames(snp.mat);
	snp.name <- rownames(snp.mat);

	###snps
	dim.snps <- dim(snp.mat);

	snps <- as.raw( as.matrix(snp.mat ) );
	snps <- array(snps, dim=dim.snps);
	colnames(snps) <- sub.name;
	rownames(snps) <- snp.name;
	class(snps) <- "SnpMatrix";

	r <- write.plink( file.base=snp.file.base, snp.major = F, snps=t(snps),
	    	id=sub.name,
	    	father=rep(0,dim.snps[2]),
	    	mother=rep(0,dim.snps[2]),
	    	sex=rep(0,dim.snps[2]),
	    	phenotype=rep(-9,dim.snps[2]),
			chromosome=chromosome,
			genetic.distance=position,
			position= position,
			allele.1 = rep("A",dim.snps[1]),
			allele.2 = rep("B",dim.snps[1]),
			na.code=0);

	cat("Genotype files have been converted into PLINK binary format(bed/bim/fam)\n");

	return(list(file.plink.bed = paste(snp.file.base, ".bed", sep=""),
   	    	file.plink.bim = paste(snp.file.base, ".bim", sep=""),
   	    	file.plink.fam = paste(snp.file.base, ".fam", sep="")));
}

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
	if (!all(ids.set == PF.gen$ind.list$member[,2]) )
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

	return(PF.gen);
}

get_snp_mat<-function( PF.gen, idx, snp.impute="mean" )
{
	gen.name <- PF.gen$gen.list$names[idx];
	snps_finding <- unique(PF.gen$gen.list$snps[which(PF.gen$gen.list$snps[,1] == gen.name), 2] );
	snp.mat <- NULL;

	if(!is.null(PF.gen$snp$matrix))
	{
		snp.mat <- get_snpstat_mat(PF.gen$snp$matrix, snps_finding, snp.impute)
	}

	if(is.null(snp.mat))
	{
		idx.range <- c(idx-50, idx+50);
		if (idx.range[1]<1) idx.range[1] <- 1
		if (idx.range[2]>PF.gen$gen.list$len) idx.range[2] <- PF.gen$gen.list$len;

		gen.names <- PF.gen$gen.list$names[idx.range[1]:idx.range[2]];
		snps <- PF.gen$gen.list$snps[which(PF.gen$gen.list$snps[,1] %in% gen.names), 2]

		PF.gen$snp$matrix<- extract_plink_data(PF.gen$options, unique(snps), PF.gen$ind.list$member);
		snp.mat <- get_snpstat_mat(PF.gen$snp$matrix, snps_finding, snp.impute)
	}

	if(!is.null(snp.mat))
		snp.mat$name <- gen.name;

	return(snp.mat);
}


read_gen_plink<-function( file.plink.bed, file.plink.bim, file.plink.fam, file.gene.set, plink.path)
{
	tb.gen <- read.table(file.gene.set, header=F, stringsAsFactors=F);
	gen.names <- unique(tb.gen[,1]);
	gen.list <- list( len=NROW(gen.names), names=gen.names, snps=tb.gen);

	tb.fam <- read.table(file.plink.fam, header=F, stringsAsFactors=F);
	ind.list <- list( member=tb.fam[,c(1,2)], removed=c() )

 	n.snp <- get_large_file_lines( file.plink.bim);

	snp <- list()
	snp$fam <- as.data.frame(tb.fam);
	options <- list( plink.path=plink.path, file.plink.bed=file.plink.bed,  file.plink.bim=file.plink.bim, file.plink.fam=file.plink.fam );

	if( n.snp * 1.0 * NROW(tb.fam) < 50*1000*2000 )
	{
		library(snpStats);
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


extract_plink_data<-function(options, snps, individuals)
{
	tmp <- tempfile( tmpdir = getwd() );
	tmp.bed <- paste(tmp, "bed", sep=".");
	tmp.bim <- paste(tmp, "bim", sep=".");
	tmp.fam <- paste(tmp, "fam", sep=".");
	tmp.all <- paste(tmp, "*", sep=".");

	file.snps <- tempfile(tmpdir = getwd() );
	file.member <- tempfile(tmpdir = getwd() );

cat("SNPS=", length(snps), "SAMPLE=", NROW(individuals), "\n")

	write.table(snps, file=file.snps, quote=F, row.names=F, col.names=F, sep="\t");
	write.table(individuals, file=file.member, quote=F, row.names=F, col.names=F, sep="\t");

	plink.cmd <- paste(options$plink.path, "--bed",  options$file.plink.bed, "--bim",  options$file.plink.bim, "--fam",  options$file.plink.fam, "--extract", file.snps, "--keep", file.member, "--make-bed", "--out", tmp, sep=" ");
cat(plink.cmd, "\n");
	plink.status <- system(plink.cmd, wait=T, intern=T);

#show(plink.status);

	library(snpStats);
	snp.mat <- snpStats::read.plink( tmp.bed,  tmp.bim, tmp.fam );

	idx.fam <- match(individuals[,2], snp.mat$fam$member)
	snp.mat$fam$member <- snp.mat$fam$member[idx.fam];
	snp.mat$genotypes <- snp.mat$genotypes[idx.fam, ];

	unlink(c(tmp.all, tmp.bed,  tmp.bim, tmp.fam, file.member, file.snps));

	return(snp.mat);
}

clone_plink_refer<-function(PF.gen)
{
	PLINK.refer <- getRefClass("PLINK.refer");
	PF.gen2 <- PLINK.refer(
				gen.list=PF.gen$gen.list,
				ind.list=PF.gen$ind.list,
				snp=list(fam=PF.gen$snp$fam),
				options=PF.gen$options);
	return(PF.gen2);
}
