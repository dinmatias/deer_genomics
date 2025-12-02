# can write the path of R to use here
# e.g., /bin/R

# NOTE
# this script is used to concatenate the different block alignments in MAF

#### USAGE
# Rscript mafConcat.R path_inputMaf path_outputMaf 
# path_inputMaf:  path to maf file to be concatenated
# path_outputMaf: path of output maf file

# Example:
# Rscript mafConcat.R /Users/dinma/Downloads/Species1to5andOG.maf C:/Users/dinma/Downloads/finalMaf.maf

# capture input
inputs <- commandArgs(trailingOnly = TRUE)

#############################################
#### parse input from command line       ####
#############################################
# parse the various inputs to different R objects

# dataDir <- "D:/tempFiles/amilFinal/"
inputPath <- as.character(inputs[1])

# outDir <- "D:/tempFiles/"
outputPath <- as.character(inputs[2])

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

blockIn <- lapply(X = 1:length(blockSt), function(x, st, en, maf){
  maf[(st[x]+1):(en[x]-1)]
}, st = blockSt, en = potEnd, maf = temp)


lenBlock <- unlist(lapply(blockIn, function(x) length(x)))

spTemp <- strsplit(x = blockIn[[which(lenBlock == 6)[1]]],
                   split = " ") 

spList <- unlist(lapply(spTemp, function(x) x[2]))

spList <- gsub(pattern = "\\.[0-9].*$", replacement = "", spList)

tempOut <- lapply(blockIn, function(x){
  
  tempList <- unlist(lapply(strsplit(x = x,
                                     split = " "), function(x) x[2]))
  
  seqList <- unlist(lapply(strsplit(x = x,
                                    split = " "), function(x) x[length(x)]))
  
  tempList <- gsub(pattern = "\\.[0-9].*$", replacement = "", tempList)
  
  # match(spList, tempList)
  seqList <- seqList[match(spList, tempList)]
  
  seqList[is.na(seqList)] <- paste0(x = rep("-", nchar(seqList[!is.na(seqList)][1])), collapse = "")
  
  tempList <- tempList[match(spList, tempList)]
  
  tempList[is.na(tempList)] <- spList[is.na(tempList)]
  
  return(list(spname = tempList, seq = seqList))
})

tempOut <- tempOut[lenBlock > 4]

seqCon <- vector(mode = "list", length = length(spList))

for(i in 1:length(tempOut)){
  
  for(sp in 1:length(spList)){
    
    if(i == 1){
      
      seqCon[[sp]] <- tempOut[[i]]$seq[sp]
      
    }else{
    
      seqCon[[sp]] <- paste0(c(seqCon[[sp]], tempOut[[i]]$seq[sp]), collapse = "")
      
    }
    
    
  }
  
  
}

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
