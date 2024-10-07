# load packages
library("easypackages") # can also install multiples with packages()
x<-c("tidyverse",  "docopt") #"reshape2", "broom", "data.table", "tibble",
libraries(x); rm(x)
options(width=220)

print("------> start of PWY TABLE COLLATION <------")
doc <- "Usage: my_program.R [options] 
 --help              show this screen
 --binDIR PATH       hard path to bin dir [default: 'wdir/TEDreadsTAX/bin']
 --suffix NAME       suffix within the file names that needs to be removed [default: '.txt']
 --outDIR NAME       path to & name of output directory to create [default: 'wdir/TEDreadsTAX/merged_tables']"
# print(doc)

args <-docopt(doc)
# print(args)
#### for manual prep #this is if the .sh does not provide arguments
# args <-docopt(doc, args = c(
#  "--binDIR" , "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/PWY_tests/test2/pwyBIN_KEGG",               # "/path/to/PWYprofiles/bin/",
#  "--suffix" , ".txt", 
#  "--outDIR" , "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/PWY_tests/test2/" ) )
args[6:8]
# print("-----------> docopt settings------^^^^")


# #=============== start scripts ==============
input=list() #instead of input=NA
input$path=args$binDIR
input$files=list.files(input$path, full.names = T, pattern = args$suffix) 
# print(paste("detecting input files from ", input$path) )
# input$files
fnames <- basename(input$files) %>% sub(args$suffix, "", .)
# input$scclist <- read.table(args$sccList, header = F, sep="\t") %>% setNames("SampNames") %>%  arrange(SampNames)
names(input$files) <- fnames;
print("------> naming input files:")
input$files

#load data from files
print("------> load data from files")
tiers=paste0(rep("Tier"), seq(1,18,1)) #assuming max(metaCyc PWY Tiers) = 18 (its actually 14)
data=list()
data=lapply(input$file, function(x){
  y<-read.table(x, header=F, check.names = F, sep=c("\t"), na.strings =c("NA",""," "), row.names = NULL,
    fill=T, quote = "", col.names = c("Counts", tiers) ) %>% remove_rownames() %>%
    group_by(across(any_of(tiers))) %>%  relocate("Counts") %>%
    unite("allTiers", any_of(tiers), sep="; ", remove=F, na.rm=T) %>%
    relocate("allTiers", .after=last_col())
})
# data[[1]]%>% tail

### grep column for MetaCyc: Tier10, for KEGG: Tier5
get_columns_with_string <- function(data, string) {
  # Get column names where any row contains the specified string
  cols_with_string <- colnames(data)[apply(data, 2, function(x) any(grepl(string, x)))]
  return(cols_with_string)
}

# head(data[[1]])
print("----> auto-detect criteria: pwyCOLUMN selected: ")
pwyCOLUMN <- get_columns_with_string(data[[1]], "PWY|map0")[[1]]; print(pwyCOLUMN)
if (pwyCOLUMN=="Tier5"){pwyNAME <- "KEGG"}else{pwyNAME <- "MetaCyc"}
print(cat("----> auto-detected criteria", pwyCOLUMN)); print(pwyNAME)

# pwyCOLUMN<- grepl("map0", data[[1]][[,1]]) ##there is no map0 or PWY in the input files. They krona
# if (length(pwyCOLUMN)>0){pwyNAME <- "KEGG"}else{pwyNAME<- "MetaCyc"}
# print(pwyCOLUMN); print(pwyNAME)


print("------> merging PWY data")
if (length(input$file)== 1){ merged<-data[[1]] %>% relocate("Counts", .after=last_col() )}else{
    suppressWarnings({ merged<-Reduce(function(x,y) merge(x,y, by=c(tiers, "allTiers"), all=T ), data) })   }
# head(merged);tail(merged)
colnames(merged)=c(tiers, "allTiers", names(data))
# head(merged); tail(merged)
merged <-  column_to_rownames(merged, pwyCOLUMN)
# row.names(merged) <- merged[[pwyCOLUMN]]##paste0("PWY", row.names(merged))
merged = merged %>%  relocate(c(any_of(tiers), "allTiers"), .after=last_col() ) %>% 
  discard(~all(is.na(.x))) %>% mutate_if(is.numeric, ~replace(. , is.na(.), 0))
# print("------> checking merged object:")
head(merged %>% select(-allTiers)); #tail(merged); dim(merged);

# output=list()
# output$dir<- args$outdir
# print(paste("------> exporting PWY data at:", output$dir))
# dir.create(output$dir)

# write.table(merged, file.path(output$dir, "merged_Counts+PWY+allTiers.txt"),
#             sep="\t", col.names=NA, quote = F, na="")

noHier <- merged %>% select(-allTiers)
write.table(noHier, file.path(args$outDIR, paste0("merged_avePWY.", pwyNAME, "_table.txt") ),
            sep="\t", col.names=NA, quote = F, na="")

Lin <- select(merged, "allTiers")
# write.table(Lin, file.path(output$dir, "merged_allTiers.txt"),
#             sep="\t", col.names=NA, quote = F, na="")

noTAX <- merged %>% select(-any_of(tiers)) %>% rename(taxonomy="allTiers")
write.table(noTAX, file.path(args$outDIR, paste0("merged_avePWY.", pwyNAME, "_4biom.txt") ),
            sep="\t", col.names=NA, quote = F, na="")
