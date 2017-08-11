library(Phoenix.datatools)
setwd("C:/Users/FiggLab4357/Desktop/Clinical Pharmacology/Kyelee Fitts/DTX Genotype Association")


# read in DMET file 
raw_data <- read.csv( "XL1_ncapk+genotype.csv", 1, stringsAsFactors = FALSE )
all_snps <- cbind ( "PatNum" = raw_data$CHP_filename , "snpID" = raw_data$dbSNP.RS.ID, "Call" = raw_data$Call , raw_data [ ,19:24 ], stringsAsFactors = FALSE )

#re-identify Patient Numbers to match Phoenix dataset's format
patnums <- all_snps$PatNum
patnums_short <- paste ( substr( patnums, 1, 1 ), substr( patnums, 7, 8 ), sep = "" )
all_snps <- cbind ( patnums_short, all_snps ) 

#read in Phoenix dataset 
dataset <- read.csv( "Docetaxel+TQTA+genotypes+sigSNPs.csv" )

# Match patnums for all_snps to patnums in dataset 
dataset_patnums <- dataset$Ptnum
dataset_patnums <- as.character(dataset_patnums)

for ( i in 1:length( dataset_patnums ) ) { 
  dataset_patnums [ i ] <- replace ( dataset_patnums [ i ], nchar( dataset_patnums [ i ] ) == 2, paste ( substr( dataset_patnums [ i ], 1, 1 ), "0", substr( dataset_patnums [ i ], 2, 2 ), sep = "" ) )
}

dataset <- cbind ( "new_patnums" = dataset_patnums, dataset )




# Add sigSNPs to dataset 

sigsnps <- snpnames

sigsnps_frame <- all_snps [ all_snps$snpID %in% sigsnps, ]

new_dataset <- new_dataset [ new_dataset$new_patnums %in% sigsnps_frame$patnums_short, ]

for ( i in 1:length( sigsnps ) ) {
  snp_i <- sigsnps_frame [ sigsnps_frame$snpID == sigsnps [ i ], ]
  snp_i$patnums_short <- as.character( snp_i$patnums_short )
  patnum_counts <- as.data.frame( table( new_dataset$new_patnums ) ) [ 25:50, ]
  patnum_counts [ , 1 ] <- as.character( patnum_counts[ , 1 ] )
  patnum_counts [ , 2 ] <- as.numeric( patnum_counts[ , 2 ] )

  snp_values <- numeric() 
  
  for ( k in 1:length( snp_i ) ) {
    for ( j in 1:length( patnum_counts ) ) {
      if ( patnum_counts [ j, 1 ] == snp_i$patnums_short [ k ] ) {
        snp_values <- append ( snp_values, rep ( snp_i$Call [ k ], patnum_counts [ j, 2 ] ) )
      }
    }
  }
  new_dataset [ , sigsnps [ i ] ] <- snp_values 
}




# Just going to do it manually. 

sigsnps_frame <- sigsnps_frame [ , c( 1, 3, 4 ) ]
sorted_snps <- data.frame( "Patnums" = unique ( sigsnps_frame [ , 1 ] ) )
for ( i in 1:length( sigsnps ) ) {
  snps_i <- sigsnps_frame [ sigsnps_frame$snpID == sigsnps [ i ], ]
  sorted_snps <- cbind ( sorted_snps, snps_i$snpID, snps_i$Call ) 
}

write.csv( sorted_snps, file = "sorted_snps.csv" )


sigsnps_data <- read.csv( "Docetaxel+genotypes+sigSNPs.csv" )

sigsnps_data$rs2276299 <- as.numeric(sigsnps_data$rs2276299)
sigsnps_data$rs2108622 <- as.numeric(sigsnps_data$rs2108622)
sigsnps_data$rs11584174 <- as.numeric(sigsnps_data$rs11584174)
sigsnps_data$rs55802895 <- as.numeric(sigsnps_data$rs55802895)
sigsnps_data$rs1051740 <- as.numeric(sigsnps_data$rs1051740)
sigsnps_data$rs13265049 <- as.numeric(sigsnps_data$rs13265049)
sigsnps_data$rs6980478 <- as.numeric(sigsnps_data$rs6980478)
sigsnps_data$rs743616 <- as.numeric(sigsnps_data$rs743616)
sigsnps_data$rs2292334 <- as.numeric(sigsnps_data$rs2292334)

write.csv( sigsnps_data, "Docetaxel+genotypes+sigsnps+numeric_covars.csv")




