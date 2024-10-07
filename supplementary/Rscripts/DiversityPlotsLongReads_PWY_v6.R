#! /usr/bin/env Rscript
### load packages
## info fo fossil/chao1: https://www.rdocumentation.org/packages/fossil/versions/0.4.0/topics/chao1
# other packages not needed for this script
#y=c("ape", "broom","readr", "reshape2", "dplyr", "tidyr", "tidyverse",
#"Biostrings","devtools","tibble", "stringr", "htmltools")
library(easypackages) #can also install multiples with packages()
x<-c("data.table", "fossil", "ggplot2", "vegan", "ampvis2" , "tidyverse", "S4Vectors")
suppressMessages({ libraries(x) }); rm(x)
options(width=220)

print("------> start of PWY Diversity Plots <------")
print("------> docopt settings ------")
library(docopt)
### prep arguments for docopt #
#(use only one space between --argument & PATH/NAME, ok to have multiple  after  or TABs)
doc <- 'Usage: DiversityPlots_v3.R [options]
  --help            show this screen.
  --mTABfile NAME      give path & name to input merged file [default: profiles/merged_tables/merged_Counts+PWY.txt ].
  --mfile NAME      give name of mapping file [default: mapping.txt ].
  --outdir NAME     give path & name of output directory to create [default: pwd/profiles/DiversityPlots].' # -> doc
doc
args <- docopt(doc); 
print(args)
## manual args set:
# args <- docopt(doc, args = c(
#    "--mTABfile"  ,  "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/PWY_tests/merged_avePWY.kegg_table.txt",
#    "--mfile"  ,  "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/PWY_tests/mapping.txt",
#    "--outdir" ,  "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/PWY_tests/DivPlots/")) #this is if the .sh does not provide arguments
# print(args)
# print("-----------> docopt settings------^^^^")

# set workdir
# setwd(args$wdir)

# #load("PWY_DiversityPlots_v3.RData")
# # 
# # #=============== starting tab input ==============
path=list()
path$outdir= paste0(args$outdir); dir.create(path$outdir); 

# #importing metadata objects
print("-----------> reading in mapping file:")
tables=NA
tables$meta <- read.table(paste0(args$mfile), header = F, sep="\t", quote = ""); #, row.names = NULL
colnames(tables$meta) <- c("SampleID", "FastqFileName", "TreatmentGroup")
row.names(tables$meta) = tables$meta$SampleID # tables$meta=tables$meta %>% relocate(SampleID)
tables$meta=tables$meta[order(tables$meta$SampleID), ]
head(tables$meta); print("-----------> done importing mapping file^^^")

### Beta diversity with Ampvis2
#https://madsalbertsen.github.io/ampvis2/articles/ampvis2.html
print("-----------> prepping AMPVIS2 object")
    # tables$amp=cbind(tables$counts, tables$tax)# or 
    # amp <- amp_load(otutable=tables$counts, metadata= tables$meta, taxonomy= tables$tax ) ## or
tables$amp=read.table(paste0(args$mTABfile),  header = T, row.names = 1, sep = "\t", quote = "" )
tables$amp$Tier0=row.names(tables$amp)
tables$amp = tables$amp %>% dplyr::rename(Species = Tier0, Genus = Tier4, Order = Tier2, Phylum = Tier1) %>%
  dplyr::select(-c(starts_with("Tier") ) ) #bc AWS::ampvis2.7.17 makes issue with the extra columns from metaCyc DB
head(tables$amp)

amp=ampvis2::amp_load(tables$amp, tables$meta)
# print("-----------> checking some AMPVIS2 data")
# head(amp$abund); tail(amp$abund)
# head(amp$tax); tail(amp$tax)
# head(amp$metadata)
# print("-----------> checking AMPVIS2 data ^^^")

print("-----------> making PCOA plot")
try({
pcoa_bray <- amp_ordinate(amp, filter_species = 0.01, type = "PCOA", sample_label_by = "SampleID", sample_label_size = 3,
                          distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity=1.5,
                          detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
#pcoa_bray$plot= pcoa_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8) )
pdf(paste0(path$outdir, "/PWY_BetaDiv_PCoA.pdf"), width=12, height=10); plot(pcoa_bray$plot); dev.off()
png(paste0(path$outdir, "/PWY_BetaDiv_PCoA.png"), width=760, height=640); plot(pcoa_bray$plot); dev.off()
})

 print("-----------> making nMDS plot")
try({
nmds_bray <- amp_ordinate(amp, filter_species = 0.01, type = "NMDS", sample_label_by = "SampleID", sample_label_size = 3,
                          distmeasure = "bray", sample_color_by = "SampleID", sample_point_size = 4, opacity = 1.5,
                          detailed_output = TRUE, transform = "none", sample_shape_by = "TreatmentGroup")
#nmds_bray$plot=nmds_bray$plot + scale_shape_manual(values=c(17, 15,13,18, 4,5,6,7,8))
pdf(paste0(path$outdir, "/PWY_BetaDiv_nMDS.pdf"), width=12, height=10); plot(nmds_bray$plot); dev.off()
})

print("-----------> making heatmap profile graph")
try({
amp_htmp=amp_heatmap(amp, group_by = "SampleID", facet_by = "TreatmentGroup", color_vector= c("darkcyan", "lawngreen"),
                     tax_aggregate = "Species", tax_show=35, plot_values = T,  normalise=T,
                     plot_values_size = 3, showRemainingTaxa = F,  tax_add = "Genus") + 
          theme(axis.text.x = element_text(angle = 45, size=10, vjust = 1),
                axis.text.y = element_text(size=8), legend.position="right")
  wratio= 8 + length(unique(amp$metadata$TreatmentGroup)) + round( (1/4)*length(amp$metadata$SampleID), 0) 
pdf(paste0(path$outdir, "/PWY_Profile_Heatmap_top35PWYs.pdf"), width=wratio, height = 10); plot(amp_htmp); dev.off()
 })


#cleanup
suppressWarnings({
rm(amp_htmp, nmds_bray, pcoa_bray)
})

#save R image
# save.image(paste0(args$outdir, "/DiversityPlots_v5.RData") )

#Attempting DiffAbund: not all dataset will have replicates & groups rich enough for statistical analysis