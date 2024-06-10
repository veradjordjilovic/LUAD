
library(tidyverse)
library(MatchIt)
library(penalizedclr)
#read in the data
load("~/Dropbox/WorkingProjects/MultiOmics/LUAD/LUAD.rdata")

# create 3 year survival variable 
# 1 for those alive after 3*365 days, 0 for dead, and NA for censored prior to 3*365 mark)
LUAD_three_year_surv <- LUAD%>% 
  mutate(three_year_surv = case_when(time > 3*365 ~ 1, 
                                     status == 0 ~ NA,
                                     status == 1 ~ 0)) 
# check the new variable
table(LUAD_three_year_surv$three_year_surv, useNA = "ifany")

temp <- cbind(LUAD$time, LUAD$status, 
              LUAD_three_year_surv$three_year_surv)
head(temp)
table(LUAD_three_year_surv$three_year_surv)

# remove NA
LUAD_no_NA <- LUAD_three_year_surv %>% filter(!is.na(three_year_surv))

# perform matching (exact on gender, and based on mahalanobis distance for age and smoking history)
res.match <- matchit(three_year_surv ~ age_clinical + 
                     tobacco_smoking_history_clinical,
                     data = LUAD_no_NA,
                     exact = "gender_MALE_clinical",
                     method = "nearest", ratio = 1,
                     distance = "mahalanobis")

# alternative matching without smoking
# res.match <- matchit(three_year_surv ~ age_clinical,
#                      data = LUAD_no_NA,
#                      exact = "gender_MALE_clinical",
#                      method = "nearest", ratio = 1,
#                      distance = "mahalanobis")
# check matching
res.match$match.matrix
LUAD_no_NA$tobacco_smoking_history_clinical[c(124, 17)]
row.names(res.match$match.matrix)

# remove unmatched subjects
LUAD_no_NA_match <- cbind(LUAD_no_NA, stratum = res.match$subclass)
LUAD_no_NA_matched <- LUAD_no_NA_match[!is.na(res.match$subclass),]

# select cnv and mrna
LUAD_select <- LUAD_no_NA_matched%>% dplyr::select(ends_with("_rna"), 
                                                   ends_with("_cnv"))
# select 1000 cnv and mrna based on association with smoking (t test)
# based on data where information on 3 year survival is not available
# data_select <- LUAD %>% 
#   mutate(nonsmoking = ifelse(tobacco_smoking_history_clinical==1, 1, 0))%>%
#   filter(time < 3*365 & status == 0)


  data_select <- LUAD %>% 
   mutate(stageIII = ifelse(stage_event_pathologic_stage_Stage.III_clinical==1, 1, 0))%>%
   filter(time < 3*365 & status == 0)

mrna_or_cnv <- LUAD %>%  
  filter(time < 3*365 & status == 0) %>%
  dplyr::select( 
  ends_with("_rna"),
  ends_with("_cnv"))
 

  
t.s <- function(x) {
  #t.res <- t.test(x ~ data_select$nonsmoking, na.action=na.omit)
  t.res <- t.test(x ~ data_select$stageIII, na.action=na.omit)
  return(t.res$p.value)
}

t.stat <- apply(mrna_or_cnv, 2, t.s) 
filtered_omics <-  names(sort(t.stat)[1:1000])
table(gsub(".*_([1-z]+)$","\\1", filtered_omics))
include_ind <- colnames(LUAD_no_NA_matched) %in% filtered_omics
LUAD_no_NA_matched_filtered_omics <- LUAD_no_NA_matched[, include_ind]



# perform the analysis
set.seed(124568)
pf  <- default.pf(response = LUAD_no_NA_matched$three_year_surv, 
                  stratum = LUAD_no_NA_matched$stratum,
                  penalized = LUAD_no_NA_matched_filtered_omics, 
                  p = table(gsub(".*_([1-z]+)$","\\1", 
                                 names(LUAD_no_NA_matched_filtered_omics))),
                  type.step1 = "comb")

lambda <- find.default.lambda(response = LUAD_no_NA_matched$three_year_surv, 
                              stratum =  LUAD_no_NA_matched$stratum,
                              penalized = LUAD_no_NA_matched_filtered_omics, 
                              p = table(gsub(".*_([1-z]+)$","\\1", 
                                             names(LUAD_no_NA_matched_filtered_omics))),
                              pf.list = list(c(4,1)))

fit.luad <- penalized.clr( response = LUAD_no_NA_matched$three_year_surv, 
                           stratum = LUAD_no_NA_matched$stratum, 
                           penalized = LUAD_no_NA_matched_filtered_omics, 
                           p = table(gsub(".*_([1-z]+)$","\\1", 
                                          names(LUAD_no_NA_matched_filtered_omics))), 
                           lambda = c(40,10),
                           standardize = TRUE,
                           alpha = 0.6)
# non-zero estimated coefficients
coeffs <- fit.luad$penalized[fit.luad$penalized!=0]


start.time <- Sys.time()
stable_res <- stable.clr.g(response = LUAD_no_NA_matched$three_year_surv,
                           stratum = LUAD_no_NA_matched$stratum, 
                           penalized = LUAD_no_NA_matched_filtered_omics, 
                           p = table(gsub(".*_([1-z]+)$","\\1", 
                                          names(LUAD_no_NA_matched_filtered_omics))),
                           standardize = TRUE,  
                           lambda.list = list(c(7,4), c(15,5), c(4,8), c(2,6)), B=50,
                           alpha = 0.6)


end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
time.taken


# Annotation
library("biomaRt")
mart <- useDataset("hsapiens_gene_ensembl", 
                   useMart("ensembl"))
ensembl_ids <- gsub("_[1-z]+$","",
                    names(LUAD_select))
df <- getBM(filters ="ensembl_gene_id", 
            attributes = c( "ensembl_gene_id" ,
                            "hgnc_symbol"), 
            values = ensembl_ids, mart = mart)

#non-zero coefficients
coeffs <- fit.luad$penalized[fit.luad$penalized != 0] 
#extracting info for plotting
library(stringi)
type <- stri_sub(names(coeffs), from = 17, to = 20)
ensebml_id <- stri_sub(names(coeffs), from = 1, to = 15)
#exract only estimated genes from the annotation
df_match <- df[df$ensembl_gene_id%in% ensebml_id,]


df_match3<- unique(df_match) 
df_match2 <- df_match3 %>% mutate(sel_prob = ifelse(ensembl_gene_id %in% 
                                                      stri_sub(selected, 
                                                               from = 1, to = 15), 1, 0))
#construct data frame for plotting
data_plot <- data.frame(coeffs = coeffs, 
                        gene = df_match2$hgnc_symbol,
                        type = type,
                        stable = df_match2$sel_prob==1)
#order
data_plot <- data_plot[order(data_plot$coeffs, decreasing = FALSE), ]
#remove genes with no hgnc match
data_plot <- data_plot[data_plot$gene != "", ]


require(ggplot2)
ggplot(data_plot, aes(x = gene, y = coeffs, col = stable))  +
  geom_point(size = 3) +
  geom_segment(aes(x = gene, xend = gene, y = 0, yend = coeffs)) +
  coord_flip() +
  #next line is to make gene order not alphabetical but according to coefficients
  scale_x_discrete(limits = unique(data_plot$gene)) +
  theme(text = element_text(size = 18))+
  scale_color_manual(values=c("red", "darkgreen"))

# Selection probability plot
plot(stable_res$Pistab, pch =20,
     col = rep(c("red", "darkgreen"), 
                                  each = table(gsub(".*_([1-z]+)$","\\1",
                                             names(LUAD_no_NA_matched_filtered_omics)))),
     ylab = "Selection probability",
     xlab = "Features")
sum(stable_res$Pistab>0.55)
abline(h = 0.55, lty = 2, pch =2) 
legend("topleft", legend=c("cnv", "mRNA"), fill=c("red", "darkgreen"), bty = "n")

selected <- names(which(stable_res$Pistab>0.55))









