import os
from os.path import join as pjoin
from os.path import exists
import pandas as pd
import csv
import snakemake.utils as sm
import yaml


############ cluster config ###############
def default_clust(cc):
    """adds __default__ values from cluster config to other
       rules lacking those keys"""
    dkeys = set(cc["__default__"].keys())
    rules = [k for k in cc if k != "__default__"]
    for r in rules:
        toadd = dkeys.difference(cc[r].keys())
        if toadd:
            cc[r].update(dict((k, cc["__default__"][k]) for k in toadd))
    return(cc)

######## load clusterFILE  ##########
with open(config["clusterFILE"]) as f:
    clust_conf = default_clust(yaml.safe_load(f))
#########  end of this mysterious section   ##### 




####### link TEMP DIR  ##################
# tmpdir = clust_conf["__default__"]["tmpdir"]  ## doesnt get read in
# print("test:", tmpdir)

JOB_ID = os.environ.get('SLURM_JOB_ID')
TMP_DIR = '/lscratch/' + JOB_ID + "/"
shell.prefix("export TMPDIR={TMP_DIR}; set -xe;")


######### read sampNames from map_file ######
### using mapFile to read in sampNames && fastqNames
### since pipeline needs map_file anyway
def readIN_mappingFile(filename):
    samples = []
    fastqNames = []
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        for line in file:
            columns = line.strip().split('\t')  # Assuming tab-separated values
            samples.append(columns[0])
            fastqNames.append(columns[1])
    return samples, fastqNames


# def read_sampList_file():
#     with open('sampList.txt') as f:
#       samples = f.read().splitlines()
#     return samples
# samples=read_sampList_file()
############# end of sampName section ################



############## environment variables ############# 
print("TMP_DIR set to:",TMP_DIR)


filename = config["map_file"]
samples, fastqNames = readIN_mappingFile(filename)
print("SampNames are:", samples)
print("fastqNames are:", fastqNames)

wDIR = os.getcwd()
print("working directory:", wDIR)



# def get_var(bash_var):
#     var = os.environ.get(bash_var)
#     return int(var)

# threads_per_job = get_var('SLURM_CPUS_PER_TASK')
# print("threads available:", threads_per_job)

# mem = os.environ.get('SLURM_MEM_PER_NODE')
# if isinstance(mem, str) and len(mem) > 0:
#     total_mem_gb = int(get_var('SLURM_MEM_PER_NODE') / 1024) 
# else:
#     total_mem_gb = int(get_var('SLURM_MEM_PER_CPU') * threads_per_job  / 1024 )
# print("memory available:", total_mem_gb, "Gb")

# tmpdir = str(get_var('SLURM_JOB_ID'))
# print("jobID:", tmpdir)




### not understanding the function
# def make_biom_file(profiles, taxfiles, samples, map_file, biom_file_name, tax_col, prof_col):
#     """
#     Make biom file from checkm profile and checkm tree_qa taxonomy across
#     all samples.

#     Args:
#         profiles (str): profile filename template with ``{sample}`` substring to format.
#         taxfiles (str): bin_TAX filename template with ``{sample}`` substring to format.
#         samples (list): list of sample names
#         map_file (str): mapping filename
#         biom_file_name (str): output biom filename
#         tax_col (str): column from bin_TAX to use
#         prof_col (str): column from profiles to use

#     Returns:
#         nothing, but produces biom file named *biom_file_name*.
#     """
#     def merge_single_profile_tax(profile, taxfile, sample, tax_col=tax_col, prof_col=prof_col):
#         """
#         Merge profile.txt and bin_TAX.txt file for single sample, and sum ``prof_col`` by ``tax_col``.
#         Returns `pandas.DataFrame` with only 2 columns with tax_col column name replaced
#         by sample name.
#         """
#         prof = pd.read_csv(profile, sep='\t')
#         tax = pd.read_csv(taxfile, sep='\t')
#         tax = tax.astype({tax_col: str})
#         tab = pd.merge(prof, tax, how="left", on="Bin Id")
#         tab = tab[[tax_col, prof_col]]
#         tab = tab.groupby(tax_col).sum()
#         tab.rename(columns={prof_col : sample}, inplace=True)
#         return(tab)

#     ## merge tables over all samples; format string to ensure files match
#     tablist = [merge_single_profile_tax(profiles.format(sample=s), taxfiles.format(sample=s), s) for s in samples]
#     if len(samples) > 1:
#         mergedtable = reduce(lambda l,r: pd.merge(l,r, how="outer", on=tax_col), tablist)
#     else:
#         mergedtable = tablist[0]
#     mergedtable = mergedtable.rename(index={'nan':'unclassified'}).fillna(0)

#     ## make biom https://biom-format.org/documentation/table_objects.html
#     observ_ids = [f"O{i}" for i in range(mergedtable.shape[0])]
#     observ_metadata = [{ 'taxonomy' : x.split(';')} for x in mergedtable.index]
#     biomtab = biom.table.Table(data=mergedtable.to_numpy(),
#                                observation_ids=observ_ids,
#                                sample_ids=list(mergedtable.columns),
#                                observation_metadata=observ_metadata,
#                                type="Taxon table")

#     ## read in mapping file
#     sample_met = biom.parse.MetadataMap.from_file(map_file)
#     biomtab.add_metadata(sample_met, axis='sample')
#     write_biom_table(biomtab, fmt="hdf5", filepath=biom_file_name)



# def clean_output():
#     misspwys = "missing_pathways.txt",
#     pwydivpl= TAX_DIVERS_PLOTS + "tmp_mapping_file.txt",
#     # asm_done=" ".join(expand(ASM + "*.done", sample=samples))
#     # qc_plots_done=" ".join(expand(QC_PLOTS + "*.done", sample=samples))
#     # pwy_plots_done=PWY_DIVERSITY_PLOTS + "*.done"
#     # tax_plots_done=TAX_DIVERSITY_PLOTS + "*.done"
#     # kraken2_clf_klogs=TED_READS_TAX_KLOGS

#     command = f"""
#     set +o pipefail
#     echo "Cleaning output..."
#     rm -rf {misspathways} 
#     """
#     return command
