##Program generate alluvial plot for the identified repurposed drugs
## Data contains experimentally validated drugs, identified repurposed drugs, and Drug type of identified repurposed drugs

#installed.packages("alluvial")
library("alluvial")


data=read.csv("alluvial_input.csv", sep=",")
getpng(filename="ncov.png", width=1100, height=500)
alluvial(data[,1:3], freq=data$Frequency, col=data$Colors, blocks=FALSE, alpha=0.8, cex=1.0, gap.width=0.4, border = "grey")
dev.off()