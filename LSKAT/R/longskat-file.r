check_plink_file<-function( file.plink.bed, file.plink.bim, file.plink.fam )
{
	cat("Checking PLINK file......\n");
	cat("* BED file =", file.plink.bed, "\n");
	cat("* BIM file =", file.plink.bim, "\n");
	cat("* FAM file =", file.plink.fam, "\n");
	
	library(snpStats);

	snp.mat <- try( read.plink( file.plink.bed,  file.plink.bim, file.plink.fam) );
	if(class(snp.mat)=="try-error")
	{
		return(list(bSuccess=F));
	}
	
	n.idv <- NROW(snp.mat$fam)
	n.snp <- NCOL(snp.mat$genotypes)
	cat("* Individuals =", n.idv, "SNPs=", n.snp, "\n")
	cat("* PLINK loading successfully.\n")
	
	return(list(bSuccess=T, family=snp.mat$fam));
}

# Phenotype File( CSV Format, Header=T  Comma seprated ) 
#
# Format 1:
# shareid, PHE1, PHE2, ...
#
# Format 3:
# shareid, TIME1, TIME2, ...

check_pheno_file<-function( file.phe.long, file.phe.time, file.plink.fam ) 
{
	cat("Checking phenotype file......\n");
	cat("* PHE.LONG =", file.phe.long , "\n");

	phe.long <- try( read.csv(file.phe.long, header=T, stringsAsFactors=F) );
	if (class(phe.long)=="try-error")
	{
		cat("! Can not open file(", file.phe.long, ")\n");
		return(list(bSuccess=F));
	}
	else{
		cat("  First 5 Items:\n");
		show(head(phe.long, n=5));
	}	
	
	cat("* PHE.TIME =", file.phe.time , "\n");
	phe.time <- NULL;
	if(!is.null(file.phe.time))
	{
		phe.time <- try( read.csv(file.phe.time, header=T, stringsAsFactors=F) );
		if (class(phe.time)=="try-error")
		{
			cat("! Can not open file(", file.phe.time, ")\n");
			return(list(bSuccess=F));
		}
		else
		{
			cat("  First 5 Items:\n");
			show(head(phe.time, n=5));
		}
	
		idx.inter <- intersect( phe.long[,1], phe.time[,1] );
		if( !( length(idx.inter)==NROW(phe.long) && length(idx.inter)==NROW(phe.time) ) )
		{
			cat("! PHE.LONG don't have consistent IDs with PHE.TIME.\n" );
			return(list(bSuccess=F));
		}		
	}
	
	cat("* PLINK.FAM =", file.plink.fam , "\n");
	phe.fam <- try( read.table(file.plink.fam, header=F, stringsAsFactors=F) );
	if (class(phe.fam)=="try-error")
	{
		cat("! Can not open file(", file.phe.fam, ")\n");
		return(list(bSuccess=F));
	}
	else{
		cat("  First 5 Items:\n");
		show(head(phe.fam, n=5));
	}
	
	m.phe <- match( as.character( phe.long[,1] ), phe.fam$V2  );
	if (length(which(is.na(m.phe))) > 0 )
	{
		cat("!", length(which(is.na(m.phe))), "IDs can not be found in the PLINK file.\n" );
		cat("! First 5 IDs:\n");
		show( head(m.phe[ which(is.na(m.phe)), ], n=5) );
	}	

	m.snp <- match( as.character( phe.fam$V2) , as.character(phe.long[,1]) );
	if (length(which(is.na(m.snp))) > 0 )
	{
		cat("!", length(which(is.na(m.snp))), " IDs can not be found in the phenotype file.\n" );
		cat("! First 5 IDs:\n");
		show( head( phe.fam[is.na(m.snp), ], n=5) );
	}	

	all.na <- which( is.na( rowSums(phe.long[,-1], na.rm=T)	) )
	if (length(all.na)>0 )
	{
		cat("!", length(all.na), "IDs dont have non-missing data.\n" );
		cat("! First 5 IDs:\n");
		show( head( phe.long[all.na, ], n=5) );
	}	

	return(list(bSuccess=T))
}

# Covariate File( CSV Format, Header=T  Comma seprated  ) 
# Format:
# 
# shareid, COV1,....

check_covariate_file<-function( file.phe.cov, file.plink.fam, y.ncov ) 
{
	cat("Checking covariate file......\n");

	cat("* COV.FILE =", file.phe.cov , "\n");
	phe.cov <- try( read.csv(file.phe.cov, header=T, stringsAsFactors=F) );
	if (class(phe.cov)=="try-error")
	{
		cat("! Can not open file(", file.phe.cov, ")\n");
		return(list(bSuccess=F));
	}
	else{
		cat("  First 5 Items:\n");
		show(head(phe.cov, n=5));
	}	
	
	if (is.na(y.ncov)) y.ncov <- NCOL(phe.cov) - 1;
	if( NCOL(phe.cov) < 1 + y.ncov)
	{
		cat("! Insufficient covariates in the covariate file, ", NCOL(phe.cov), "<", 1 + y.ncov, ".\n" );
		return(list(bSuccess=F));
	}	

	phe.fam <- try( read.table(file.plink.fam, header=F, stringsAsFactors=F) );
	if (class(phe.fam)=="try-error")
	{
		cat("! Can not open file(", file.phe.fam, ")\n");
		return(list(bSuccess=F));
	}
	
	m.phe <- match( as.character(phe.cov[,1]), as.character(phe.fam$V2)  );
	if (length(which(is.na(m.phe))) > 0 )
	{
		cat("!", length(which(is.na(m.phe))), "IDs, can not be found in the PLINK file.\n" );
		cat("! First 5 IDs:\n");
		show( head(phe.cov[is.na(m.phe), ], n=5 ) );
	}
	
	m.snp <- match( as.character(phe.fam$V2) , as.character(phe.cov[,1]) );
	if (length(which(is.na(m.snp))) > 0 )
	{
		cat("!", length(which(is.na(m.snp))), "IDs, can not be found in the phenotype file.\n" );
		cat("! First 5 IDs:\n");
		show( head(phe.fam[is.na(m.snp), ], n=5 ) );
	}
		
	all.na <- which( is.na( rowSums(phe.cov[,-1], na.rm=T)	) )
	if (length(all.na)>0 )
	{
		cat("!", length(all.na), "IDs dont have non-missing data.\n" );
		cat("! First 5 IDs:\n");
		show( head( phe.cov[all.na, ], n=5) );
	}	

	return(list(bSuccess=T))
}

# GENE SET File( Header=F, 2 columns, Space seprated ) 
# Format:
# 
# snpname, gene.name

check_geneset_file<-function( file.gene.set ) 
{
	cat("Checking gene definition file......\n");
	cat("* GEN.SET.FILE =", file.gene.set , "\n");

	tb <- try( read.table(file.gene.set, sep=" ", header=F) );
	if(class(tb)=="try-error")
		return(list(bSuccess=F));

	genes <- unique(tb[,1]);
	
	cat("* GENEs =", length(genes), "\n");
	cat("* SNPs =", NROW(tb), "\n");
	
	return(list(bSuccess=T))
}
