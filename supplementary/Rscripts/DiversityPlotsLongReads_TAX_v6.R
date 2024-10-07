#! /usr/bin/env Rscript
### load packages
## info fo fossil/chao1: https://www.rdocumentation.org/packages/fossil/versions/0.4.0/topics/chao1
library(easypackages) #can also install multiples with packages()
x<-c("data.table", "fossil", "ggplot2", "vegan", "ampvis2" , "tidyverse")
# other packages not needed for this script
#y=c("ape", "broom","readr", "reshape2", "dplyr", "tidyr", "tidyverse",
#"Biostrings","devtools","tibble", "stringr", "htmltools")
suppressMessages({ libraries(x) }); rm(x)
options(width=220)

print("-----------> docopt settings---------------------")
library(docopt)
### prep arguments for docopt #
#(use only one space between --argument & PATH/NAME, ok to have multiple  after  or TABs)
doc <- 'Usage: DiversityPlots_v1.R [options]
  --help            show this screen.
  --mTABfile NAME      give path to input file [default: merged_asmbTAX_table_wLin.txt ].
  --mfile NAME      give name of mapping file [default: mapping.txt ].
  --outdir NAME     give path & name of output directory to create [default: pwd/profiles/DiversityPlots].' # -> doc
doc
args <- docopt(doc); 
print(args)
## manual args set:
# args <- docopt(doc, args = c(
#    "--wdir"   ,  "/Users/angelovaag/Documents/Nephele_pipelines/WmGS_v2/tests",
#    "--mTABfile"  ,  "TAXprofiles/merged_asmbTAX_table_wLin.txt",
#    "--mfile"  ,  "mapping.txt",
#    "--outdir" ,  "TAXprofiles/DiversityPlots")) #this is if the .sh does not provide arguments
# print(args)
# print("-----------> docopt settings------^^^^-----------")

# #load("TAX_DiversityPlots_v1.RData")
# # 
# # #=============== starting tab input ==============
path=list()
path$outdir= paste0(args$outdir); dir.create(path$outdir); 

# #------------------------------- alpha and beta diversity -------------------------------------
# #importing counts, meta and tax objects
tables=NA

print("-----------> reading in MAPPING file")
tables$meta <- read.table(paste0(args$mfile), header = T, sep="\t", quote = "", comment=""); #, row.names = NULL
names(tables$meta)[1] <- "SampleID"
# colnames(tables$meta) <- tables$meta[1, ] ##c("SampleID", "FastqFileName", "TreatmentGroup")
row.names(tables$meta)= tables$meta$SampleID # tables$meta=tables$meta %>% relocate(SampleID)
tables$meta=tables$meta[order(tables$meta$SampleID), ]
tables$meta$TreatmentGroup <- as.factor(tables$meta$TreatmentGroup)
Ncols<- ncol(tables$meta); Nrows <- nrow(tables$meta)
try({ head(tables$meta, 10) })

print("-------------> read in COUNTS-TAX table (LongReads version)")
tables$merged <- read.table(paste0(args$mTABfile), header = T, row.names = 1, sep="\t", na.strings = "Unassigned", fill=T) 
try({ tables$merged[1:Ncols, 1:Nrows]  }) #tables$merged[1:Ncols, 1:Nrows] 
linCol <- "Lineage"
taxRanks <- c("kingdom","phylum", "class", "order", "family",  "genus", "species")

print("-------------> extracting COUNTS (LongReads version)")
tables$counts <- tables$merged %>% select(-all_of(c(linCol, taxRanks) ) ) %>% as.data.frame()
try({ print(tables$counts[1:5, 1:5])  })

print("-------------> extracting TAX (LongReads version)")
tables$tax <- tables$merged %>% select(all_of(taxRanks))
try({ print(tables$tax[8:15, 3:7]) }) 

print("-----------> Calculating alpha diversity indexes")
library(vegan)
alpha1 <- matrix(NA, ncol = 6, nrow = Nrows) %>% as.data.frame()
row.names(alpha1) <- tables$meta$SampleID 
colnames(alpha1) <- c("nReads","SpRichness", "Evenness", "Shannon", "InvSimpson", "chao1")
alpha1$nReads =  colSums(tables$counts)
alpha1$SpRichness=vegan::specnumber(tables$counts, MARGIN = 2)
alpha1$Shannon=vegan::diversity(as.matrix(tables$counts), "shannon", MARGIN = 2)
alpha1$InvSimpson=vegan::diversity(tables$counts, "inv", MARGIN = 2)
alpha1$Evenness=alpha1$Shannon/log(specnumber(tables$counts, MARGIN = 2))
alpha1$chao1=apply(tables$counts, 2, chao1) #library(fossil)
try({ print(alpha1[1:6, 1:Nrows]) }) 


print("-----------> Adding alpha indices to MAPPING");
# alpha1$temp=cbind(alpha1$SpRichness, alpha1$Evenness, alpha1$reads, alpha1$Shannon, alpha1$InvSimpson, alpha1$chao1)#, tree_div)
# colnames(alpha1$temp)=c("SpRichness", "Evenness", "NumbReads", "Shannon", "InvSimpson", "chao1")
# tables$meta=cbind(tables$meta, round(alpha1$temp, 2)); #tables$meta
tables$meta <- cbind(tables$meta, round(alpha1, 2)); 

print("-------------> meta with alphas (LongReads version)")
try({ head(tables$meta, 10) }) 
write.table(tables$meta, paste0(path$outdir, "/TAX_AlphaDiv.txt"), sep="\t", col.names=NA, quote=F, na="")


#Grouped boxplots for alpha diversity
sst=NA
sst$meta=tables$meta
sst$meta$group1=tables$meta$SampleID
sst$meta$group2=tables$meta$TreatmentGroup
sst$meta=sst$meta[, c("group1","group2",  "SpRichness","Evenness" , "Shannon", "InvSimpson")]
sst$melt=melt(as.data.table(sst$meta), id.var=c("group1", "group2")); #(sst$melt)
sst$melt$group1=factor(sst$melt$group1)

print("-----------> plotting alpha diversity")
try({
  p1 <- ggplot(sst$melt, aes(group2, value)) + 
    geom_boxplot(aes(x=group2, y=as.numeric(as.character(value), stat="identity",
                                            position="dodge"), fill=group2)) + theme_bw() +
    facet_wrap(~variable, scales="free") + labs(y="Alpha Diversity Measures", x="TreatmentGroup", fill="TreatmentGroup") +
    theme(text=element_text(size=16), axis.text.x = element_text(angle=45, hjust=1, size=10),
          legend.position = "right", legend.title=element_text(size=12),
          legend.text=element_text(size=12))
  pdf(paste0(path$outdir , "/TAX_AlphaDiv", ".pdf"), width=12 ) ; print(p1); dev.off()
})
sst=NULL;


# ### Beta diversity with Ampvis2
# #https://madsalbertsen.github.io/ampvis2/articles/ampvis2.html
print("-----------> prepping AMPVIS2 object - longReads")
# tables$amp=cbind(tables$counts, tables$tax)# or 
amp <- amp_load(otutable=tables$counts, metadata= tables$meta, taxonomy= tables$tax ) ## or
# tables$amp=read.table(paste0(path$tables, "/merged_Counts+TAX.txt"),  header = T, row.names = 1, sep = "\t", quote = "" )
# try({tables$amp[1:5,1:5] })
# amp=amp_load(tables$amp, tables$meta)
# print("-----------> checking some AMPVIS2 data")
# try({  amp$abund[1:5, 1:5] })
# try({    amp$tax[1:5, 3:7] })
# try({amp$metadata[1:3,1:3] })
# print("-----------> checking AMPVIS2 data")

print("-----------> making PCOA plot")
try({
  title=paste("PCoA ordination plot of TAX content, using Bray-Curtis distances")
  pcoa_bray <- amp_ordinate(amp, filter_species = 0.01, type = "PCoA", sample_label_by = "SampleID", sample_label_size = 3,
                            distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity=1.5,
                            detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
  #pcoa_bray$plot= pcoa_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8) ) + labs(title=title)
  pdf(paste0(path$outdir, "/TAX_BetaDiv_PCoA.pdf"), width=12, height=10); print(pcoa_bray$plot); dev.off()
  png(paste0(path$outdir, "/TAX_BetaDiv_PCoA.png"), width=760, height=640); print(pcoa_bray$plot); dev.off()
})

print("-----------> making nMDS plot")
try({
  title=paste("nMDS ordination plot of TAX content, using Bray-Curtis distances")
  nmds_bray <- amp_ordinate(amp, filter_species = 0.01, type = "NMDS",  sample_label_by = "SampleID", sample_label_size = 3,
                            distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity = 1.5,
                            detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
  #nmds_bray$plot=nmds_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8))  + labs(title=title)
  pdf(paste0(path$outdir, "/TAX_BetaDiv_nMDS.pdf"), width=12, height=10); plot(nmds_bray$plot); dev.off()
})

print("-----------> making heatmap profile graph")
try({
  amp_htmp=amp_heatmap(amp, group_by = "SampleID", facet_by = "TreatmentGroup",
                       tax_aggregate = "Species", tax_show=35, plot_values = T,  normalise=T,
                       plot_values_size = 3, showRemainingTaxa = T, tax_add = "Class") +
    theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
          axis.text.y = element_text(size=8), legend.position="right")
  wratio= 8 + length(unique(amp$metadata$TreatmentGroup)) + round( (1/4)*length(amp$metadata$SampleID), 0) 
  pdf(paste0(path$outdir, "/TAX_Profile_Heatmap.pdf") , width=wratio, height = 10 ); plot(amp_htmp); dev.off() ##
}) ## made plot size flexible

print("-----------> making rarecurve graph")
try({
  # print(colSums(amp$abund) %>% sort() )
  rmin <- min(unlist(lapply( colSums(amp$abund), \(x)  x[x>1000] ) ) ) ##  x[x!=0] ## selects first value >1K reads from nREADs/sample
  xmin=2*rmin #x-axis min-limit at 2x 1K+ readDepth value (or at least 2*lowest-non-0 sample)
  xmax=max(sort(colSums(amp$abund)))
  rareAT <- round(xmax/2, 0) ## rarifying to shorten calcTime
  suppressMessages({ amp2 <- amp_rarefy(amp, rarefy = rareAT) }) ## saving on calcTime
  step_size <-round(rareAT/100,0) ## round(xmax/100, 0) #
  suppressWarnings({
    rrc <- amp_rarecurve(amp2, color_by = "SampleID", facet_by = "TreatmentGroup") +#, stepsize = step_size) +  #facet_by = "TreatmentGroup", 
      xlim(0,xmin) + theme(legend.position="bottom") + # legend.text = element_text(size=6) ) +  
      guides(color=guide_legend(ncol=10) ) +  #ncol=6, nrow=10
      labs(title=paste("Rarecurve"), x="Number reads (seqDepth)", y="Number observed ASVs") 
    # plot(rrc)
  }) 
  hratio= 10+ length(unique(amp$metadata$TreatmentGroup))
  wratio= 10+ length(unique(amp$metadata$TreatmentGroup))
  pdf(paste0(path$outdir, "/TAX_RareficaitonCurve.pdf"), width=wratio, height=hratio); plot(rrc); dev.off() #, width=12, height=10, height=hratio
  
}) # made rarecurves more flexible in terms of how they are run and plot size

print("-----------> making rank abundance / Whittaker plot")
try({
  rab=amp_rankabundance(amp, group_by = "TreatmentGroup", log10_x = T)
  pdf(paste0(path$outdir, "/TAX_RankAbundanceCurve.pdf")); plot(rab); dev.off()
})

#cleanup
suppressWarnings({
  rm(amp_htmp, amp2, nmds_bray, p1, pcoa_bray, sst, rrc, rab)
})

#save R image
#save.image(paste0(args$outdir, "/DiversityPlots_v3.RData") )

#Attempting DiffAbund
#Need replicates