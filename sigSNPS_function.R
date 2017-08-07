library( ggplot2 )
library(reshape2)

sig.snps <- function ( directory, significance, exportplots = FALSE ) {  
  
  setwd( directory )
  
  raw_data <- read.csv( "XL1_ncapk+genotype1.csv", 1, stringsAsFactors = FALSE )
  
  snp.frame <- cbind ( "PatNum" = raw_data$CHP_filename , "snpID" = raw_data$dbSNP.RS.ID, "Call" = raw_data$Call , raw_data [ ,19:24 ], stringsAsFactors = FALSE )
  
  unique.snps<- unique ( snp.frame$snpID )
  
  
  # Function to obtain p-value for ANOVA
  
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  sigsnps1 <- data.frame ()
  sigsnps2 <- data.frame ()
  
  for ( j in 4:9 ) { 
    for ( i in 1:length ( unique.snps ) ) { 
      test.frame <- snp.frame[ snp.frame$snpID == unique.snps [ i ], ]
      
      # When there are more than 3 genotypes, use ANOVA 
      
      if ( length ( unique ( test.frame$Call ) ) >= 3 ) { 
        anova <- lm ( test.frame[ , j ]~Call, data = test.frame )
        pvalue <- lmp ( anova )
        if ( pvalue <= significance ) {
          print ( paste ( "SNP = ", unique.snps[ i ], ":: Parameter = ", colnames ( test.frame ) [ j ], ":: p = ", pvalue ) ) 
          plot.frame <- cbind ( "Call" = test.frame$Call, "Parameter" = test.frame[ , j ] )
          plot.frame <- data.frame (plot.frame)
          plot.frame[ , "Call2"] <- cbind ( as.character( plot.frame [ , 1 ] ) )
          plot.frame[ , "Parameter2" ] <- cbind ( as.numeric( levels( plot.frame [ , 2 ] ) )[ plot.frame [ , 2 ] ] )
          title <- test.frame$snpID [ 1 ]
          ylabel <- colnames ( test.frame ) [ j ] 
          print ( ggplot ( plot.frame, aes ( x = Call2, y = Parameter2 ) ) + geom_point() + geom_boxplot( ) + ggtitle ( title ) + ylab ( ylabel ) )
          
          if ( exportplots ) {
            ggsave( paste( title, ".jpeg", sep = "" ) )
          }
          
          sigsnps1 <- rbind ( sigsnps1, cbind( title, colnames ( test.frame ) [ j ] ) ) 
        }
        
      # When there are exactly 2 genotypes with at least 2 data values per genotype, use t-test. 
        
      } else if (  length ( unique ( test.frame$Call ) ) == 2 & sum ( test.frame$Call == unique ( test.frame$Call ) [ 1 ] ) > 2 & sum ( test.frame$Call == unique ( test.frame$Call ) [ 2 ] ) > 2  ) { 
          ttest <- t.test ( test.frame[ , j ]~test.frame$Call )
          if ( ttest$p.value <= significance ) {
            print ( paste ( "SNP = ", unique.snps[ i ], ":: Parameter = ", test.frame [ , j ], ":: p = ", ttest$p.value ) ) 
            plot.frame <- cbind ( "Call" = test.frame$Call, "Parameter" = test.frame[ , j ] )
            plot.frame <- data.frame (plot.frame)
            plot.frame[ , "Call2"] <- cbind ( as.character( plot.frame [ , 1 ] ) )
            plot.frame[ , "Parameter2" ] <- cbind ( as.numeric( levels( plot.frame [ , 2 ] ) )[ plot.frame [ , 2 ] ] )
            title <- test.frame$snpID [ 1 ]
            ylabel <- colnames ( test.frame ) [ j ] 
            print ( ggplot ( plot.frame, aes ( x = Call2, y = Parameter2 ) ) + geom_point() + geom_boxplot( ) + ggtitle ( title ) + ylab ( ylabel ) )
            
            if ( exportplots ) {
              ggsave( paste( title, ".jpeg", sep = "" ) )
            }
            
            sigsnps2 <- rbind ( sigsnps2, cbind( title, colnames ( test.frame ) [ j ] ) )       
        }
      } 
    }
  }
  
  # Collecting significant SNPs in a data frame
  
  sigsnps <- rbind( sigsnps1, sigsnps2 )
  colnames( sigsnps ) <- c ( "SNPs", "Parameter" )
  return ( sigsnps )
}
      
list.snps <- sig.snps ( directory = "~/Desktop/Clinical Pharmacology/Kyelee Fitts/Docetaxel_phoenix", 
                        significance = 5 * 10^-5, 
                        export = TRUE )





snpnames <- list.snps [ , 1 ]
snpnames <- as.character ( unique(snpnames) )
sort(snpnames)
snpnames

