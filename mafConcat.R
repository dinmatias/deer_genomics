# can write the path of R to use here
# e.g., /bin/R

# NOTE
# this script is used to concatenate the different block alignments in MAF

#### USAGE
# Rscript mafConcat.R path_inputMaf path_outputMaf 
# path_inputMaf:  path to maf file to be concatenated
# numPresent:     the minimum number of lineages aligned to include block
# path_outputMaf: path of output maf file

# Example:
# Rscript mafConcat.R /Users/dinma/Downloads/Species1to5andOG.maf 5 C:/Users/dinma/Downloads/finalMaf.maf

# capture input
inputs <- commandArgs(trailingOnly = TRUE)

#############################################
#### parse input from command line       ####
#############################################
# parse the various inputs to different R objects

# dataDir <- "D:/tempFiles/amilFinal/"
inputPath <- as.character(inputs[1])

numPresent <- as.integer(inputs[2])

# outDir <- "D:/tempFiles/"
outputPath <- as.character(inputs[3])



#### end inputs

# read the maf file 
# use readLines to read each line as an element of a vector
temp <- readLines(inputPath)

# find the start of block
# block starts with "a score="
# use grep to search for this pattern
blockSt <- grep(pattern = "a score=",
                x = temp)

potEnd <- grep(pattern = "^$",
               x = temp)

# remove indices prior to the first block
# that is prior to the first appearance of "a score="
potEnd <- potEnd[!(potEnd < blockSt[1])]

# pull out the different alignment blocks
blockIn <- lapply(X = 1:length(blockSt), function(x, st, en, maf){
 
  maf[(st[x]+1):(en[x]-1)]
 
}, st = blockSt, en = potEnd, maf = temp)

# determine how many lineages were aligned for each block
lenBlock <- unlist(lapply(blockIn, function(x) length(x)))

# pull out alignment block that is complete
# max lenBlock should be complete
spTemp <- strsplit(x = blockIn[[which(lenBlock == max(lenBlock)[1])[1]]],
                   split = " ") 

# pull out the species name from the complete block
# the point is to use this pattern as template
spList <- unlist(lapply(spTemp, function(x) x[2]))

# remove the additional numbers added by last
spList <- gsub(pattern = "\\.[0-9].*$", 
               replacement = "", spList)

# use the spList as template
# arrange the alignments to correspond to the template
# add a gap for lineages with no alignment
tempOut <- lapply(blockIn, function(x){
  
  tempList <- unlist(lapply(strsplit(x = x,
                                     split = " "), 
                            function(x) x[2]))
  
  seqList <- unlist(lapply(strsplit(x = x,
                                    split = " "), 
                           function(x) x[length(x)]))
  
  tempList <- gsub(pattern = "\\.[0-9].*$", replacement = "", tempList)
  
  # match(spList, tempList)
  seqList <- seqList[match(spList, tempList)]
  
  seqList[is.na(seqList)] <- paste0(x = rep("-", nchar(seqList[!is.na(seqList)][1])), collapse = "")
  
  tempList <- tempList[match(spList, tempList)]
  
  tempList[is.na(tempList)] <- spList[is.na(tempList)]
  
  return(list(spname = tempList, seq = seqList))
})

# remove blocks with less than numPresent aligned lineages
tempOut <- tempOut[lenBlock > numPresent]

# holder for the concatentaed sequence
# blocks are being concatenated by linage
seqCon <- vector(mode = "list", length = length(spList))

# concatenate
for(i in 1:length(tempOut)){
  
  for(sp in 1:length(spList)){
    
    if(i == 1){
      
      seqCon[[sp]] <- tempOut[[i]]$seq[sp]
      
    }else{
    
      seqCon[[sp]] <- paste0(c(seqCon[[sp]], tempOut[[i]]$seq[sp]), collapse = "")
      
    }
  }
}

# write the output
for(i in 1:length(seqCon)){
  
  if(i == 1){
    
    write("##maf version=1 scoring=multiz", outputPath, append = T)
    write("", outputPath, append = T)
    write("a score=620137250.0", outputPath, append = T)
    
  }
  write(paste0(c("s", " ", spList[i], " ", "255699", " ", 
                 "676528", " ", "+" , " ", "676528", " ", 
                 seqCon[[i]]), collapse = ""), 
        outputPath, append = T)  
}
