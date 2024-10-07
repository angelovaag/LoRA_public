# load packages
library("easypackages") # can also install multiples with packages()
x<-c("tidyverse","docopt") #"reshape2",  #"data.table", "broom",  "tibble", , "stats"
suppressMessages(libraries(x)); rm(x)
options(width=220)

#prep arguments for docopts
doc <- "Usage: my_program.R [options]
--help           show this screen
--binDIR PATH    hard path to bin dir [default: 'getwd()/TEDreadsTAX/bin']
--suffix NAME    the suffix to be removed from the file names to obtain only the SampName
--outTAB NAME    path & name of output directory to create [default: 'getwd()/TEDreadsTAX/merged_tables']"
# doc
args <- docopt(doc)
# args
## #for manual prep
#  args <-docopt(doc, args = c(
#   "--binDIR" , "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/TAX_tests/test_taxmerge/kronabin1/",
#   "--suffix"    , "_4krona.txt",
#   "--outTAB"     , "/Users/angelovaag/Documents/Nephele_pipelines/minION_longRead/Rscript_tests/TAX_tests/test_taxmerge/merged_asmbTAX_table.txt" ))
# args


# #=============== starting script ==============
input=list()
input$files=list.files(args$binDIR, full.names = T) #pattern = args$infix,
# print("------> detecting input files:")
# input$files
fnames <- basename(input$files) %>% sub(args$suffix, "", .)
names(input$files) <- fnames;

# print("------> naming input files:")
# input$files

#load data from files
print("------> load data from files")
data=list()
krColNames <- c("contigName", "Counts", "TAXid")

data=lapply(input$files, function(x){
  y<-read.table(x, header=F, check.names = F, sep=c("\t"), na.strings =c("NA",""," "),
                fill=T, quote = "", col.names = krColNames ) %>% select(-krColNames[1]) %>%
    group_by(TAXid) %>% reframe(sumCounts=sum(Counts))
    return(y)
})
print("------> checking data loading objects:")
head(data[[1]]);# tail(data[[1]])
# names(data)

print("------> merging data") ### will work with only 1 file as well
suppressWarnings({ merged=Reduce(function(x,y) merge(x,y, by="TAXid", all=T ), data) }  )
colnames(merged) <- c("TAXid", fnames)
merged <- merged %>% mutate_if(is.numeric, ~replace(. , is.na(.), 0)) # will not be making rownames 
  # column_to_rownames("TAXid") 
# head(merged)


print("------> checking merged object:")
head(merged)

print(paste("------> exporting table to", args$outTAB))

write.table(merged, file.path(args$outTAB),
            sep="\t", col.names=T, quote = F, na="", row.names = F)
