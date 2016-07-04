#   Original author:
#   Liana Nice, Saint Paul, MN
impute <- function (ScaffoldData,
                    ParentData,
                    Maps,
                    WhichMap = "M11",
                    SNP_names="iSelect", 
                    Contig = NA,
		                Variant = NULL,
                    fams = 1:length(Families), 
                    AltMap = NA) {
  
  ##ScaffoldData needs columns: Line, Family, Markers, Markers should correspond to SNP_names
  #### All data should be relative to common parent.
  
  ##ParentData needs rows = markers corresponding to SNP_names OR Contig, and col names = parent names and SNP_names or Contig.  
  #### Common parent should not be included, and all data should be relative to common parent.
  
  ##Maps should have a column of names == "SNP_names", and two columns: Chrom_"WhichMap" and cM_"WhichMap"
  
  ##AltMap is if you want a secondary ordering based on two maps
  
  ############ Running Code:   #####################################################
  ### Establishing maps and sorting by maps 
  if(is.na(Contig)) Contig <- SNP_names
  MapChrom <- paste("Chrom_", WhichMap, sep="")
  Map_cM <- paste("cM_", WhichMap, sep="")
  Maps[,MapChrom]<-as.factor(Maps[,MapChrom])
  Maps <- droplevels(Maps)
  Chroms <- levels(Maps[,MapChrom])
  Parents <- colnames(ParentData)[!colnames(ParentData) %in% c(Contig, Variant)]
  ScaffoldData <- droplevels(ScaffoldData)
  Families <- levels(ScaffoldData[["Family"]])
  
  NumFams <- length(fams)
  NumChrom <- length(Chroms)
  
  ## Order of data by map
  if(length(Families) != length(Parents)) {
    ErrorFrame <<- cbind (Families, Parents) 
    stop("Number of families and parents are not equivalent.  View object ErrorFrame.") }
    
  	if(!is.na(AltMap)) {
  	  AltMapChrom <- paste("Chrom_", AltMap, sep="")
  	  AltMap_cM <- paste("cM_", AltMap, sep="")
  	  MapOrder <-order(Maps[,MapChrom], Maps[,Map_cM],Maps[,AltMapChrom], Maps[,AltMap_cM], Maps[,Contig])
  	} else {
  	  MapOrder <-order(Maps[,MapChrom], Maps[,Map_cM], Maps[,Contig])  
  	}   
  
  Maps <- Maps[MapOrder,]
  Maps <- Maps[!is.na(Maps[,Map_cM]),]
  
  MapOrder <- match(names(ScaffoldData)[-1:-2], Maps[,SNP_names])
  ScaffoldData <- ScaffoldData[,c(TRUE,TRUE,!is.na(MapOrder))]
  MapOrder <- na.omit(MapOrder)
  MapOrder <- order(MapOrder)
  ScaffoldData <- ScaffoldData[,c(1:2,(MapOrder+2))]

  MapOrder <-  match(ParentData[[Contig]], Maps[[Contig]])
  ParentData <- ParentData[!is.na(MapOrder),]
  MapOrder <- na.omit(MapOrder)
  MapOrder <- order(MapOrder)
  ParentData <- ParentData[MapOrder,]
  if(!is.null(Variant)) { rownames(ParentData) <- ParentData[[Variant]]
  contig_variant <- data.frame(ParentData[[Contig]], ParentData[[Variant]])
  } else {rownames(ParentData) <- ParentData[[Contig]]}
  
  TotMarkImpute <- nrow(ParentData)
  message (paste("Imputing", TotMarkImpute, "Markers"))
  ScafNames <- colnames(ScaffoldData)
  
  

  
  #  empty list of matricies 
  imputations <- vector("list", length=NumFams)
  SegMarkers <- vector("list", length=NumFams)
  MarkerFrames <- vector("list", length=NumFams)
  SegMarkerMaps <- vector("list", length=NumFams)
  
  ####  Imputation Loop ####
  
  #########    Start Family Loop 
  for (f in fams) {
  #### Establish Family and Parent data  = seg markers, scaffold markers, and parent values  
    Family <- Families[f]
    Parent <- Parents[f]
    
    message("Imputing ", Family, paste(", ", f, " of ", NumFams, ", ", Sys.time(), sep=""))
    
    ParentRow <- ScaffoldData[ScaffoldData$Line %in% Parent,]
    
    SegMarkers[[f]]  <- ScafNames[!is.na(ParentRow) & ParentRow == 1]
    SegMarkerMaps[[f]] <- Maps[Maps[,SNP_names] %in% SegMarkers[[f]],]
    MarkerFrames[[f]] <- ScaffoldData[ScaffoldData$Family == Family, c("Line", "Family",SegMarkers[[f]])]
    NumIndv <- nrow(MarkerFrames[[f]])   
    
    #### Make imputations Matrix 
    imputations[[f]] <- matrix(nrow = NumIndv, ncol= TotMarkImpute)
    colnames(imputations[[f]]) <- rownames(ParentData)
    rownames(imputations[[f]]) <- MarkerFrames[[f]]$Line
    
    
    #### Establish frames for Chromosomes within Families 
    MarkerFrameChr <- vector("list", length=7)
    SegMarkerMapsChr <- vector("list", length=7)  
    
  #### Set imputations which are not segregating in the parent = 0, subsest everything else by this  
    imputations[[f]][,ParentData[[Parent]] == 0] <- 0
    isSeg <-!is.na(ParentData[[Parent]]) & ParentData[[Parent]] != 0
    ImputeMarks <- ParentData[[Contig]][isSeg]
    
    NumImputeMarks <- length(ImputeMarks)  
    ImputeMarksCM <- match(ImputeMarks, Maps[,Contig])
    ImputeMarksCM <- Maps[ImputeMarksCM,Map_cM]

    message(" Number of individuals:", NumIndv, "\n Number of segregating markers to impute:", NumImputeMarks, 
		"\n Number of non-segregating markers:", sum(!isSeg))

    
  #### Start Chromosome Loop ####
  c = 1  
  
  for (h in 1:NumChrom) {
      Chromosome <- Chroms[h]
      message( "Chrom ", Chromosome, paste(", ", h, " of ", NumChrom, ", ", Sys.time(), sep=""))
  
      #### Establish data for Chromosomes within Families 
      SegMarkersChr <- as.character(SegMarkerMaps[[f]][SegMarkerMaps[[f]][,MapChrom] == Chromosome, SNP_names])
      MarkerFrameChr[[h]] <- MarkerFrames[[f]][,SegMarkersChr]
      row.names(MarkerFrameChr[[h]])<- MarkerFrames[[f]][ , "Line"]
      SegMarkerMapsChr[[h]] <- SegMarkerMaps[[f]][SegMarkerMaps[[f]][ , SNP_names] %in% SegMarkersChr,]
      NumSegMarkers <- length(SegMarkersChr)
      
      for (m in 1:(NumSegMarkers-1)) {
  
  #### Begining of chromosomes ####
        while (c <= NumImputeMarks & ImputeMarksCM[c] <= SegMarkerMapsChr[[h]][1,Map_cM]) {
            m2 <- MarkerFrameChr[[h]][,m]
            while (sum(is.na(m2))>0) {m2[is.na(m2)] <- MarkerFrameChr[[h]][is.na(m2),m+1]}
            imputations[[f]][,ParentData[[Contig]] == ImputeMarks[c] & isSeg] <- m2
  
          c=c+1
        }
  
  
  #### Between scaffold markers ####   
        
        while (ImputeMarksCM[c] >= SegMarkerMapsChr[[h]][m,Map_cM] & 
                 ImputeMarksCM[c] <= SegMarkerMapsChr[[h]][m+1,Map_cM] &
                   c <= NumImputeMarks) {      
          ### Set Vector of flanking marker 1 locations and calls
          LocFM1 <- rep(SegMarkerMapsChr[[h]][m,Map_cM], times=NumIndv)
          MCall1 <- MarkerFrameChr[[h]][,m]
          m1 <- m
          while (sum(is.na(MCall1))>0) {
           
            if(m1-1 < 1) {
              ErrorFrame <<- MarkerFrameChr[[h]]
              stop("No flanking marker call to to calculate interval for marker ", ImputeMarks[c],
                   " on Chrom ", Chromosome, ". \n \n", "See scaffold marker ", SegMarkersChr[m1], " in individuals: \n ", 
                   paste(rownames(imputations[[f]])[is.na(MCall1)],collapse=" ,") )
              
            }else{
            
            MCall1[is.na(MCall1)] <- MarkerFrameChr[[h]][is.na(MCall1),m1-1]
            LocFM1[is.na(MCall1)] <- SegMarkerMapsChr[[h]][m1-1,Map_cM]
            m1 <- m1-1  
          }
          }        
          
          ### Set Vector of flanking marker 2 locations and calls
          LocFM2 <- rep(SegMarkerMapsChr[[h]][m+1,Map_cM], times=NumIndv)
          MCall2 <- MarkerFrameChr[[h]][,m+1]
          m2 <- m+1
          
          while (sum(is.na(MCall2))> 0) {
            
            if(m2+1 > NumSegMarkers) {
              ErrorFrame <<- MarkerFrameChr[[h]]
              stop("No flanking marker call to to calculate interval for marker ", ImputeMarks[c],
                   " on Chrom ", Chromosome, ". \n \n", "See scaffold marker ", SegMarkersChr[m2], " in individuals: \n ", 
                   paste(rownames(imputations[[f]])[is.na(MCall2)],collapse=" ,") )
            }else{
              
              MCall2[is.na(MCall2)] <- MarkerFrameChr[[h]][is.na(MCall2),m2+1]
              LocFM2[is.na(MCall2)] <- SegMarkerMapsChr[[h]][m2+1,Map_cM]
              m2 <- m2+1
            }
          }
            
            GapTotal <- LocFM2-LocFM1 
            Gap1 <- ImputeMarksCM[c]-LocFM1
            Gap2 <- LocFM2 - ImputeMarksCM[c]
            
            #### Filling imputed values
            
          imputations[[f]][,ParentData[[Contig]] == ImputeMarks[c] & isSeg] <- (MCall1*(1-(Gap1/GapTotal)))+(MCall2*(1-(Gap2/GapTotal)))
        
          c = c+1      
        }
        
  #### Ends of Chromosomes ####  
   
        while (c <= NumImputeMarks & ImputeMarksCM[c] >= SegMarkerMapsChr[[h]][NumSegMarkers,Map_cM])
        {
          m2 <- MarkerFrameChr[[h]][,m]
          while (sum(is.na(m2))>0) {m2[is.na(m2)] <- MarkerFrameChr[[h]][is.na(m2),m-1]}
          imputations[[f]][,ParentData[[Contig]] == ImputeMarks[c] & isSeg] <- m2
          
          c=c+1
        } 
      }
    }
}
  
   imputations
}

write.imputations <- function(imputations, csv_file, list_file, Maps, SNP_names) {
  save(imputations, file=list_file)
  
  all_imputations <- do.call("rbind",imputations)
  all_imputations <- t(all_imputations)
  all_imputations <- as.data.frame(all_imputations)
  all_imputations[,SNP_names] <- row.names(all_imputations)
  all_imputations <- merge(all_imputations,Maps, by.x=SNP_names, by.y=SNP_names )
  write.csv(all_imputations, file=csv_file)

}
