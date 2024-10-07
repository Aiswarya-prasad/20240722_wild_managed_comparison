# Load required libraries

if(!require(tidyverse)){
  install.packages(pkgs = 'tidyverse', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyverse)
}

if(!require(gplots)){
  install.packages(pkgs = 'gplots', repos = 'https://stat.ethz.ch/CRAN/')
  library(gplots)
}

if(!require(RColorBrewer)){
  install.packages(pkgs = 'RColorBrewer', repos = 'https://stat.ethz.ch/CRAN/')
  library(RColorBrewer)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Missing arguments", call.=FALSE)
} else {
  input <- args[1]
  output.heatmap <- args[2]
  output.table <- args[3] 
  dist <- args[4]
}

###### Heatmap ######

### Abundance ###

# Read matrix

matabjac <- as.matrix(read.csv(input, sep = ";", header = T, row.names = 1))
colnames(matabjac) <- rownames(matabjac)
matabjac[lower.tri(matabjac)] <- t(matabjac)[lower.tri(matabjac)] #symmetrize matrix

# Transform 0-1 distances in 0-100 similarity measure
matabjac = (1 - matabjac) * 100

## Computing mini-maxi for colour palette

mini=min(matabjac[])
maxi=max(matabjac[row(matabjac)!=col(matabjac)]) # ignoring the diagonal
trueMax=max(matabjac[]) # typically the value in the diagonal = 100
q25=quantile(matabjac[row(matabjac)!=col(matabjac)],0.25,1)
q50=quantile(matabjac[row(matabjac)!=col(matabjac)],0.5,1)
q75=quantile(matabjac[row(matabjac)!=col(matabjac)],0.75,1)

n=100 # number of steps between 2 colors

## We use the quantiles to ignore some outlier values in the matrix (values<mini will have colour of mini and values>maxi will have a colour between brown and grey23)
mini=max(q25-1.5*(q75-q25),0)
maxi=min(q75+1.5*(q75-q25),trueMax)

palette=colorRampPalette(c("green", "yellow", "red", "brown", "grey23"))(n = 5*n-1)

## Checking if maxi = trueMax
trueMax.needed=ifelse(maxi<trueMax,"T","F")

if(trueMax.needed){
  breaks=c(seq(mini,maxi,length=4*n),seq(maxi+1e-5,trueMax,length=n))
  # breaks are equally distributed in the range mini-maxi (intervals can be different in the range maxi-trueMax, containing very few points)
} else {
  breaks=c(seq(mini,maxi,length=5*n))
}

# Dendrogram is obtained with the symetric matrix
distance    = dist(matabjac)
cluster     = hclust(distance, method="average")
dendrogram  = as.dendrogram(cluster)

png(output.heatmap, width = 8000, height = 8000)

heatmap.2(matabjac,
          trace = "none",
          dendrogram = "row",
          key = F,
          Rowv = dendrogram,
          Colv = rev(dendrogram),
          col = palette,
          breaks = breaks,
          labRow = rownames(matabjac),
          cexRow = 2,
          labCol = rownames(matabjac),
          cexCol = 2,
          cellwidth = 60, 
          cellheight = 60,
          margins = c(20,20),
          main = paste0(dist, " distance of kmer abundance similarity"))

# Adding the colour scale
par(fig=c(0.05,0.3,0.9,0.95), mar=rep(2,4), new=TRUE)

if(trueMax.needed){
  
  diff=maxi-mini
  breaksToMaxi=breaks[1:(4*n)] # using only breaks from mini to maxi
  black.width=max(diff/9)
  black.space=max(diff/9)
  
  plot(c(mini,maxi+black.width+black.space),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxt="n",xaxs = "i", yaxs = "i")
  rect(breaksToMaxi[-length(breaksToMaxi)],0,breaksToMaxi[-1],2,col=palette,border=NA)
  
  
  ti=pretty(breaksToMaxi)
  ti=ti[ti<maxi]
  axis(1,at=c(ti,maxi+black.space+black.width/2), label=c(ti,trueMax), cex.axis=10)
  
  # Here plotting the TrueMax colour with a white space
  rect(maxi+black.space,0,maxi+black.space+black.width,2,col=palette[5*n-1],border=NA)
  rect(maxi,-0.1,maxi+black.space,2.1,col="white",border=NA)
  
} else{
  plot(range(breaks),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxs = "i", yaxs = "i")
  rect(breaks[-length(breaks)],0,breaks[-1],2,col=palette,border=NA)
}

dev.off()

###### Top 50 most similar samples for each sample ######

### Abundance ###

subset_top <- function(matrix, sample){
  sub <- as.data.frame(matabjac[ ,sample]) # extract column corresponding to sample
  sub <- cbind(sub, rownames(sub)) # add column with sample names
  colnames(sub) <- c("similarity", "samples") # rename columns
  sub <- sub %>% arrange(desc(similarity)) # arrange by highest to lowest similarity
  return(sub$samples[1:50])
}

top <- as.data.frame(matrix(0, nrow = 50, ncol =1))
  
for (S in colnames(matabjac)){
  top <- cbind(top, subset_top(matabjac, S))
}

# Remove first empty column

top <- top[ ,-1]

# Put target sample as colname 

colnames(top) <- top[1, ]

top <- top %>%
  pivot_longer(cols = colnames(top), names_to = "assembly", values_to = "reads") %>%
  arrange(assembly)

# Write output

write.table(top, file = output.table, sep = "\t", quote = F, row.names = F, col.names = T)