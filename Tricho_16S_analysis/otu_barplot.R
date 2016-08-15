library(phyloseq)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#Load OTU data from biom file
otu_table_file <- 'denovo_otus99/otu_table_filt0001.biom'
otus <- import_biom(otu_table_file)
df <- as.data.frame(otu_table(otus))

#Convert read counts to percentages
df.norm <- sweep(df, 2, colSums(df), FUN='/') * 100
df.norm["otu"] <- row.names(df.norm)

#Convert table from 'wide' format to 'long' format for plotting
df.norm.melt <- melt(df.norm,
                     id.vars = 'otu',
                     variable.name = 'Sample',
                     value.name = 'percent')

#Make plot and save as .png file
plot_file <- 'denovo_otus99/otu_table_filt0001_barplot.png'
png(plot_file)
title <- "OTUs clustered at 99% identity"
p <- ggplot(df.norm.melt, aes(x=Sample,y=percent,fill=otu)) + geom_bar(stat='identity', col='black') + ggtitle(title) + scale_color_brewer(type = 'qual', palette = 'Set3')
print(p)
dev.off()

