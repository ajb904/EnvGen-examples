library(phyloseq)
library(ggplot2)
library(reshape2)

### Loading read data
#load our representative set of sequenced reads
otu_reads <- read.delim('denovo_otus99/otu_table_filt0001.txt', skip = 1, row.names = 1, header=TRUE)
otu_reads <- cbind(otu_reads, ref=rep(0, nrow(otu_reads)))

#Load tree
tr <- read.tree(file = 'align_otus_w_known_seqs/aligned_seqs/all_seqs_aligned.tre')

# Load taxonomy assignments for the sequenced reads, and coerce into correct format
tax.reads <- read.delim('denovo_otus99/uclust_assigned_taxonomy/Tricho_seqs_rep_set_tax_assignments.txt', header=F)

tax.reads <- with(tax.reads, cbind(V1, colsplit(V2, '; __', names=c('Kingdom',
                                                                  'Phylum',
                                                                  'Subsection',
                                                                  'Family',
                                                                  'Genus',
                                                                  'Species',
                                                                  'Strain'))))
rownames(tax.reads) <- tax.reads$V1
tax.reads <- tax.reads[,2:8]

#Set Strain = OTU name, so that we can display OTU names and known strains on the same plot
tax.reads$Strain <- row.names(tax.reads)


###Loading reference data
#Load reference taxonomy assignments, and coerce into correct format
tax.ref.raw <- read.delim('known_tricho_seqs/known_tricho_rep_set_taxSilva.txt', header=F)

tax.ref <- with(tax.ref.raw, cbind(V1, colsplit(V2, '; __', names=c('Kingdom',
                                                              'Phylum',
                                                              'Subsection',
                                                              'Family',
                                                              'Genus',
                                                              'Species',
                                                              'Strain'))))
rownames(tax.ref) <- tax.ref$V1
tax.ref <- tax.ref[,2:8]

#Do some annoying faff to turn sort out the lower levels of the SILVA taxonomy
tax.strings <- tax.ref.raw$V2
genus <- sapply(tax.strings, function(x) rev(unlist(strsplit(as.character(x), "; __")))[2])
minclass <- sapply(tax.strings, function(x) rev(unlist(strsplit(as.character(x), "; __")))[1])
pat <- "(.*?_.*?)_.*"
species <- sub(pat, "\\1", minclass)

tax.ref$Genus <- genus
tax.ref$Species <- species
tax.ref$Strain <- minclass

# Make a dummy OTU table for reference sequences
otu.ref <- matrix(rep(c(0,0,2),nrow(tax.ref)), ncol=3, byrow = T)
row.names(otu.ref) <- row.names(tax.ref)
colnames(otu.ref) <- c('Tn004', 'Tn019','ref')


###Combine read and reference data
#Combine the read and reference OTU tables and taxonomy tables
otu.comb <- otu_table(rbind(otu_reads, otu.ref), taxa_are_rows = T)
tax.comb <- tax_table(as.matrix(rbind(tax.reads, tax.ref)))

#Create a Phyloseq object from the OTU table, taxonomy table and tree
full_otu <- phyloseq(otu_table(otu.comb), tax.comb, phy_tree(tr))


###Plotting
#Plot a phylogenetic tree and save it to a pdf file
pdf('uclust_SILVA_treeplot.pdf', paper='a4r')
p <- plot_tree(full_otu, color = "Sample", size="Abundance", label.tips = "Strain", sizebase = 10, base.spacing = 0.05)
print(p)
dev.off()
