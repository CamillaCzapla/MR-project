##Two Sample MR code

##import data build37 from GWAS catalog (in terminal)
scp .\33462485-GCST90016970-EFO_0007874-Build37.f.tsv.gz camilla@10.22.9.205:/home/camilla/TwoSampleMR
##rest of code in RStudio
##making a data frame for exposure data
##possible template for later
##necessary packages for MR
library(dplyr)
library(data.table)
library(TwoSampleMR)
##read file
dt = fread("33462485-GCST90016943-EFO_0007874-Build37.f.tsv.gz")
##select rows from the read in data table
exposure = select(dt, variant_id, beta, standard_error, effect_allele, other_allele, p_value) #, effect_allele_frequency)
write.table(exposure, file = "Family_Oxalobacteraceae_exp_dat.txt", quote = FALSE, row.names = FALSE)
##sorting by p value
object = fread("Family_Oxalobacteraceae_exp_dat.txt")
attempt = arrange(object, p_value)
##slicing data table
sliced = slice_head(attempt, n=3000)
##filtering by p-value
filtered_Family_Oxalobacteraceae_exp_dat = filter(attempt, p_value < 1e-05)
write.table(filtered_Family_Oxalobacteraceae_exp_dat, file = "filtered_Family_Oxalobacteraceae_exp_dat.txt", quote = FALSE, row.names = FALSE)

#how i made the large connected table
fwrite(res, file = "Genus_Bifidobacterium_res.txt")
fwrite(res, file = "attempt_class_Actinobacteria_res_copy.txt", append = TRUE)

##THIS WORKED!!!!!!!!!!
test_data = microbiome_exp_dat %>% filter(SNP == "rs3768998" | SNP == "rs55842567" | SNP == "rs6660520")

##exposure data
microbiome_exp_dat <- read_exposure_data(
  filename = "filtered_class_Actinobacteria_exp_dat.txt",
  sep = " ",
  snp_col = "variant_id",
  beta_col = "beta",
  se_col = "standard_error",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  ## eaf_col = "effect_allele_frequency",
  pval_col = "p_value",
  ## units_col = "Units",
  ## gene_col = "Gene",
  ## samplesize_col = "n"
)

##phenotype default is "exposure"...changing it to "Gut microbiome"
microbiome_exp_dat$exposure <- "class Actinobacteria"

##LD clumping
microbiome_exp_dat <- clump_data(microbiome_exp_dat)

##outcome
ao <- available_outcomes()

##note: ulcerative colitis is 'ieu-a-970' and rheumatoid arthritis is 'ieu-a-833'

##examples of extracting SNPs
chd_out_dat <- extract_outcome_data(
  snps = microbiome_exp_dat$SNP,
  outcomes = 'ieu-a-9'
)

uc_out_dat <- extract_outcome_data(
  snps = microbiome_exp_dat$SNP,
  outcomes = 'ieu-a-970'
)

ra_out_dat <- extract_outcome_data(
  snps = microbiome_exp_dat$SNP,
  outcomes = 'ieu-a-833'
)

##HARMONISE DATA (CREATES NEW COMBINED DATA FRAME WITH EXPOSURE AND OUTCOME DATA)
dat <- harmonise_data(microbiome_exp_dat, chd_out_dat)
## for replicating ulcerative colitis and Bifidobacterium
dat <- harmonise_data(microbiome_exp_dat, uc_out_dat)
dat <- harmonise_data(microbiome_exp_dat, ra_out_dat)
##Perform MR
res <- mr(dat)
##Perform MR with specified methods
res <- mr(dat, method_list=c("mr_egger_regression", "mr_ivw"))
res <- mr(dat, method_list=c("mr_ivw"))
res1 <- mr(dat, method_list =c("mr_wald_ratio"))

dat$exposure <- "Family Oxalobacteraceae"
dat$outcome <- "Ulcerative Colitis"
res$exposure <-"Family Oxalobacteraceae"
res$outcome <- "Ulcerative Colitis"

##save data and results
fwrite(dat, file = "MR_data_Ulcerative_Colitis_vs_Family_Oxalobacteraceae.txt", quote = FALSE, sep = " ")

##Sensitivity analyses
het<-mr_heterogeneity(dat)
plt<-mr_pleiotropy_test(dat)
sin<-mr_singlesnp(dat)

##Combine results
all_res<-combine_all_mrresults(res,het,plt,sin,ao_slc=T,Exp=T,split.exposure=F,split.outcome=T)
head(all_res[,c("Method","outcome","exposure","nsnp","b","se","pval","intercept","intercept_se","intercept_pval","Q","Q_df","Q_pval","consortium","ncase","ncontrol","pmid","population")])

##Subset on method
subset_on_method(res)

subset_on_method(
  mr_res,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted"
)

##Convert log odds ratios into odds ratios with 95% confidence intervals
generate_odds_ratios(res)

##VISUALS

#Scatter plot
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
p1[[2]]

fun_plot <- mr_scatter_plot(res1, dat)
fun_plot [[1]]

##saving scatter plot as .png file
png(file="/home/camilla/TwoSampleMR/Rheumatoid_Arthritis_vs_Family_Oxalobacteraceae_IVW_scatterplot.png",
    width=600, height=350)
p1[[1]]
dev.off()

##also for png (or pdf - replace .png with .pdf)
library(ggplot2)
ggsave(p1[[1]], file="IVW_plot_Ulcerative_Colitis_vs_class_Actinobacteria_trial.png", width=9, height=7.5, units="cm")

graph = p1[[1]] 
graph_bigger_text = graph + theme_bw(20)

##the file I will save results to after replicating MiBiogen
fwrite(res, file = "MiBioGen_replicated_results.csv", append = TRUE)

#Forest plot
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
p2[[1]]


##"A report can be generated that performs all MR analyses, sensitivity analyses, and plots,
##and presents them in a single self-contained html web page, word document, or pdf document."
mr_report(dat)

##to experiment with later
single_method = subset_on_method(res)
silly_plot <- mr_scatter_plot(single_method, dat)
silly_plot[[1]]
