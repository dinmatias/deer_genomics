# can write the path of R to use here
# e.g., /bin/R

# NOTE
# this script examines concordance between estimated 

#### USAGE
# Rscript concordEstim.R output1 output2 
# output1:  path to first estimate
# output2:  path to second estimate
# output3:  a path with jpg extension for the plot

# Example:
# Rscript concordEstim.R /Users/dinma/Downloads/out.txt /Users/dinma/Downloads/out2.txt /Users/dinma/Downloads/out.jpg

# capture input
inputs <- commandArgs(trailingOnly = TRUE)

#############################################
#### parse input from command line       ####
#############################################
# parse the various inputs to different R objects

inputPath1 <- as.character(inputs[1])

inputPath2 <- as.character(inputs[2])

inputPath3 <- as.character(inputs[3])

#####

temp1 <- readLines(inputPath1)

temp1 <- temp1[grep(pattern = "^Posterior mean ", x = temp1):length(temp1)]

temp1 <- gsub(pattern = "\\s+", replacement = " ", temp1)

out1 <- strsplit(x = temp1[grep(pattern = "t_n", x= temp1)],
                 split = " ")

useTab1 <- do.call("rbind", out1)[ , 1:2]

rm(temp1)

temp2 <- readLines(inputPath2)

temp2 <- temp2[grep(pattern = "^Posterior mean ", x = temp2):length(temp2)]

temp2 <- gsub(pattern = "\\s+", replacement = " ", temp2)

out2 <- strsplit(x = temp2[grep(pattern = "t_n", x= temp2)],
                 split = " ")

useTab2 <- do.call("rbind", out2)[ , 1:2]

useTab2 <- useTab2[match(useTab1[, 1], useTab2[ ,1]), ]

useTab <- data.frame(nodes = useTab1[, 1], 
                     estim1 = useTab1[ , 2],
                     estim2 = useTab2[ , 2])

lmOut <- lm(estim1 ~ estim2, data = useTab)

lmSum <- summary(lmOut)

jpeg(inputPath3, width = 400, height = 400)
plot(x = useTab$estim2, 
     y = useTab$estim1, 
     xlab = "estim2", ylab = "estim1")

lines(x = lmOut$model$estim2,
      y = lmOut$fitted.values)

legend("topleft", legend = paste("R2 =", lmSum$r.squared),
       bty = "n")

dev.off()
