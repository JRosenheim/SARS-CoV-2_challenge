# SARS-CoV-2 Challenge Study RNAseq viral biomarkers: Data preparation and viral biomarker scores for blood mRNAseq data
# Started 2023-01-18
# Joshua Rosenheim & Rishi Gupta

library(tidyverse)
library(readxl)
library(stringr)
library(data.table)
library(pROC)
library(ggpubr)
library(rms)
library(ggcorrplot)
library(gtsummary)
library(flextable)
library(ggfortify)
library(janitor)
library(wesanderson)


# Load deduplicated TPM matrix----

tpm <- read.csv("CHIC-blood_tpm_PC0.001_log2_genesymbol_dedup.csv", row.names = 1, header = T)


# Transpose and format TPM matrix----

tpm_t <- data.frame(t(tpm)) %>% 
  rownames_to_column("sample_id") %>%
  mutate(sample_id = str_remove_all(sample_id, "_S[:digit:][:digit:][:digit:]"), # remove suffixes
         sample_id = str_remove_all(sample_id, "_S[:digit:][:digit:]"),
         sample_id = str_remove_all(sample_id, "_S[:digit:]"),
         sample_id = str_remove_all(sample_id, "X")) # remove 'X' prefix

saveRDS(tpm_t, "tpm_t.rds")
tpm_t <- readRDS("tpm_t.rds")


# Metadata----

metadata <- read.csv("CHIC-blood_metadata_no.seropositives.csv") %>%
  clean_names() %>%
  rename(sample_id = 1) %>%
  rename(subject_id=participant_id) %>%
  mutate(study_day = str_remove(study_day, "d"),
         study_day = case_when(
           study_day == "0.pre" ~ "0", 
           study_day == "0.post" ~ "0.5", 
           T ~ study_day),
         study_day = as.numeric(study_day)) %>%
  select(-c("filename", "fastq_id", "m_counts_mapped", "read_depth_strata", "m_total_reads", "percentage_pseudoaligned", "percentage_pseudoaligned_strata", "age_at_screening", "serology", "library_prep_batch", "library_prep_method", "sequencing_batch"))


# Join and check data has joined correctly----

tpm_metadata <- metadata %>%
  left_join(tpm_t) %>%
  rowwise() %>%
  mutate(female_score = 
           mean(c(
             RPS4X	,
             AMMECR1	,
             CD40LG	,
             ZRSR2	,
             EFHC2	,
             CA5B	,
             ZFX	,
             EIF1AX	,
             MORC4	,
             UBA1	,
             SYAP1	,
             DDX3X	,
             FUNDC1	,
             NKRF	,
             ZC4H2	,
             PIM2	,
             SHROOM4	,
             USP9X	,
             SMC1A	,
             NUP62CL	,
             ERCC6L	,
             NAA10), na.rm=T) -
           mean(c(KDM5D	,
                  RPS4Y1	,
                  EIF1AY	,
                  USP9Y	,
                  DDX3Y	,
                  UTY	,
                  PRKY	,
                  ZFY	,
                  CYBB	,
                  TMSB4Y), na.rm=T)) %>%
  ungroup()

tpm_metadata %>%
  ggplot() +
  geom_boxplot(aes(x=sex, y=female_score), outlier.alpha = 0) +
  geom_jitter(aes(x=sex, y=female_score), height=0.25, width=0.25)


# Save prepared data----

saveRDS(tpm_metadata, "tpm_metadata.rds")
tpm_metadata <- readRDS("tpm_metadata.rds")


# Read in signature gene list and check if in RNAseq data----

gene_list <- read.csv("Viral RNA Biomarkers - Gene list_2023-03-13.csv") %>%
  mutate(signature=str_remove_all(signature, "-"),
         gene_symbol_new=str_replace(gene_symbol_new, "-", "."),
         gene_symbol_new = case_when(
           gene_symbol_original=="KIAA1324" ~ "ELAPOR1",
           T ~ gene_symbol_new))

gene_list$present <- gene_list$gene_symbol_new %in% colnames(tpm_metadata)
gene_list$gene_symbol_new[gene_list$present==F] # Print missing genes
# [1] "Unable to find DBET"    "LOC100288289 withdrawn" "LOC100291626 withdrawn" "LOC652301 withdrawn"    "rRNA?"                  "Not found"              "lncRNA?"
gene_list_present <- gene_list %>% dplyr::filter(present==T) #Drop missing genes (7 - one Tsalik probe maps to 4 of these (DBET, LOC100288289, LOC100291626, LOC652301); 3 Steinbrink19 genes missing)


# Import modules----

STAT1_module <- read.csv("CRM_STAT1.csv", header=F) %>%
  pull(1)


# Check how many genes of each module are available----

table(STAT1_module %in% colnames(tpm_metadata)) # All 41 present


# Keep metadata and genes of interest----

tpm_sel <- tpm_metadata %>%
  dplyr::select(
    any_of(c(colnames(metadata),
             gene_list_present$gene_symbol_new,
             "CENPE", "TYMS", "MCEMP1", "CD4",
             
             # TB genes
             "BATF2",
             "ANKRD22",
             "OSBPL10",
             "GBP5",
             "DUSP3",
             "KLF2",
             
             # Modules
             STAT1_module
    )))


# Generate existing signature scores----

tpm_sel <- tpm_sel %>%
  
  ungroup() %>%
  
  rowwise() %>%
  
  mutate(Herberg2 = FAM89A - IFI44L,
         
         Pennisi2 = ADGRE1 - IFI44L,
         
         Sampson4 = (ISG15 + OASL) - (IL16 + ADGRE5),
         
         Sampson10 = (DIAPH2+GBP2+TLR5) - (IL7R+GIMAP4+FGL2) - ((ISG15+OASL) - (IL16+ADGRE5)),
         
         Sweeney7 = ((IFI27+JUP+LAX1)/3) - (((HK3+TNIP1+GPAA1+CTSB)/4) * 4/3),
         
         TrouilletAssant6	= median(c(IFI27, IFI44L, IFIT1, ISG15, RSAD2, SIGLEC1)),
         
         AndresTerre11 = mean(c(CD38, HERC5, HERC6, IFI6, IFIH1, LGALS3BP, LY6E, MX1, PARP12, RTP4, ZBP1)),
         
         Yu3 = mean(c(IL1RN, SIGLEC14, MCTP1)),
         
         Lopez7 = 
           -4*LY6E +
           -4*SIGLEC1 +
           -1*IFIT1 +
           -1*TRDV3 +
           1*ARG1 +
           4*CD177 +
           5*VNN1,
         
         Lydon15 = 
           (CADM1 *	0.048847731) +
           (CFAP45	* -0.057551661) +
           (GIT2	* -0.609978511) +
           (GLIPR1 *	-0.011838467) +
           (IFI27	* 0.063542302) +
           (IFNGR2	* -0.127951678) +
           (ELAPOR1	* -0.004273478) +
           (LAPTM4B	* -0.23099519) +
           (MX1*	0.009412439) +
           (RABGAP1L*	0.142094167) +
           (RSAD2*	0.013684041) +
           (RUNX1*	-0.148085233) +
           (SERPING1*	0.000115328) +
           (SIGLEC1*	0.314917697) +
           (SMPD1*	-0.132317095),
         
         Zaas48 = 
           IFI27	*	0.365764	+
           SIGLEC1	*	0.256661	+
           IFI44L	*	0.24436	+
           RSAD2	*	0.176117	+
           IFI44	*	0.132044	+
           ISG15	*	0.013072	+
           SPATS2L	*	0.011503	+
           LAMP3	*	0.011203	+
           SERPING1	*	0.010568	+
           IFIT1	*	0.007308	+
           OAS3	*	0.005626	+
           OASL	*	0.004477	+
           MX1	*	0.004222	+
           OAS1	*	0.003973	+
           LY6E	*	0.003591	+
           SEPTIN4	*	0.003591	+
           HERC5	*	0.00331	+
           IFIT3	*	0.003251	+
           OAS2	*	0.002795	+
           XAF1	*	0.002441	+
           ATF3	*	0.002359	+
           RTP4	*	0.002341	+
           IFIT2	*	0.002288	+
           TNFAIP6	*	0.001787	+
           GBP1	*	0.001726	+
           IFI6	*	0.001685	+
           CTSL	*	0.001654	+
           IFIT5	*	0.001332	+
           CCL2	*	0.001132	+
           DDX58	*	0.001127	+
           mean(c(LILRB1,LILRB2))	*	0.000817	+
           ADAR	*	0.000485	+
           SOCS1	*	0.000476	+
           RUBCNL	*	0.000456	+
           CUZD1	*	0.000403	+
           PRSS21	*	0.000396	+
           SOCS5	*	0.000382	+
           NOD2	*	0.000326	+
           RPL30	*	0.000313	+
           GM2A	*	0.000293	+
           HLA.DOB	*	0.00028	+
           NLRP3	*	0.000235	+
           GAPDH	*	0.000231	+
           IL16	*	0.000176	+
           ENOSF1	*	0.000092	+
           NDUFA10	*	0.000045	+
           PPIA	*	0.000037	+
           SOCS2	*	0.000017,
         
         Tsalik33 = 
           RTCB	*	0.038998	+
           LY6E	*	0.166043	+
           TNFSF10	*	0.005084	+
           TES	*	0.07874	+
           SP100	*	0.02937	+
           IRF2	*	0.074576	+
           CTBP1	*	-0.01392	+
           IRF9	*	0.034534	+
           CAMK1	*	0.111394	+
           DDX3Y	*	-0.06712	+
           CIB2	*	0.223868	+
           mean(c(DEFA1,DEFA1B,DEFA3)) 	*	-0.08786	+
           SORBS1	*	0.243737	+
           CD160	*	0.118889	+
           ARPC3	*	0.582264	+
           KPNB1	*	0.742946	+
           POLR2F	*	-0.03119	+
           OASL	*	0.185097	+
           KIDINS220	*	-0.01023	+
           RPS21	*	-0.5768	+
           BTF3	*	0.103261	+
           DAPK2	*	-0.06503	+
           SCAPER	*	0.326241	+
           ADGRD2	*	-0.00075	+
           DUX4 *-0.0343	+
           MRPS18B	*	1.07798	+
           TMEM165	*	-0.61377	+
           PPCDC	*	0.221446	+
           TRMT13	*	-0.15077	+
           ANKRD11	*	-0.40545	+
           GIMAP6	*	0.25509	+
           CFAP45	*	-0.03456	+
           ZNF335	*	-0.11226,
         
         Sweeney11 = mean(c(CEACAM1, ZDHHC19, NMRK1, GNA15, BATF, C3AR1)) - (mean(c(FAM214A,TGFBI,MTCH1,RPGRIP1,HLA.DPB1)) *6/5),
         
         Henrickson16 = mean(c(PIMREG, HESX1, PLSCR1, OASL, OAS3, RSAD2, USP18, IFI27)) - mean(c(EIF4B, EIF3L, FBL, QARS1, EEF2, UXT, NT5E, ESYT1)),
         
         Sweeney3 = ((GBP5+DUSP3)/2)-KLF2,
         
         Suliman2 = ANKRD22-OSBPL10,
         
         Buturovic6 = mean(c(DEFA4, BATF, HK3)) - mean(c(TGFBI,LY86, HLA.DPB1)), # this is an approximation as weightings not publicly available
         
         IFI27.CENPE = IFI27 + CENPE,
         
         IFI27_MX1 = (IFI27 + MX1)/2,
         
         IFI27_IFI44L = (IFI27 + IFI44L)/2,
         
         Cappuccio11 = mean(c(GUCD1, PIF1, BANF1, EHD3, TCEAL3)) - mean(c(ARAP2, SLC25A46, SLK, ROCK2, DOCK5, TVP23B)),
         
         Rao8 = mean(c(SMARCD3, ICAM1, EBI3)) - mean(c(JUP, SUCLG2, IFI27, FCER1A, HESX1)),
         
         Rao8_bacterial = mean(c(SMARCD3, ICAM1, EBI3)),
         
         Rao8_viral = mean(c(JUP, SUCLG2, IFI27, FCER1A, HESX1)),
         
         Li3 = NAGK *	0.076193749 +
           HERC6 *	-0.734076039 +
           IGF1R *	0.629199527,
         
         Steinbrink19 = 
           #MT-RNR2	*	-0.5201	+
           VPS29	*	-0.3985	+
           MMD	*	-0.1855	+
           IZUMO4	*	-0.182	+
           #AC015912.3	*	-0.1795	+
           ATP5MK	*	-0.0969	+
           TMEM170B	*	-0.0669	+
           #SNHG8	*	-0.0008	+
           CCDC71	*	0.027	+
           BTBD9	*	0.0543	+
           PBDC1	*	0.0712	+
           CMPK2	*	0.1287	+
           TMEM199	*	0.1691	+
           ISG15	*	0.2129	+
           HERC6	*	0.2211	+
           DDA1	*	0.232	+
           LY6E	*	0.5983	+
           MAGED2	*	0.603	+
           PIGT	*	0.8054,
         
         Xu2 = IFI44L	*	1.657	+
           PI3	*	-1.468,
         
         Ravichandran10 = (mean(c(HK3, GYG1, MMP9)) - mean(c(DNMT1, PRF1, MX1, IFI27, IFI44, ISG15, EPSTI1)))*3/7, #3/7 is irrelevant to ranking but included by authors
         
         Ravichandran10_gm = exp(mean(log(c(HK3, GYG1, MMP9)))) - exp(mean(log(c(DNMT1, PRF1, MX1, IFI27, IFI44, ISG15, EPSTI1))))*3/7, #Authors seem to use geometric mean despite being on log scale
         
         GomezCarballa3 = BATF	*	-1.16	+
           ISG15	*	0.64	+
           DNMT1	*	1.23
  ) %>% 
  ungroup() %>%
  mutate(
    STAT1_module = rowMeans(dplyr::select(., any_of(STAT1_module)))
  )


# Create vector of existing signatures----

gene_list_present$signature[gene_list_present$signature == "Zaas47"] <- "Zaas48" # Zaas48 incorrectly named Zaas47 in list of signatures
existing_signatures <- c(levels(factor(gene_list_present$signature)))
existing_signatures
existing_signatures <- existing_signatures[!existing_signatures %in% c("OLFM4", "Buturovic6", "Sweeney11")]
existing_signatures


# Transform all signatures to Z scores based on day 0 samples----

get_zscore<- function(variable, dataset=tpm_sel){
  x_all <- dataset[[variable]]
  x_controls <- dataset %>% filter(study_day==0) %>% pull(variable)
  mean_controls <- mean(x_controls)
  sd_controls <- sd(x_controls)
  x_zscore <- (x_all-mean_controls)/sd_controls
  dataset[[variable]] <- x_zscore
  return(dataset)
}

for(i in c(existing_signatures, "CENPE", "IFI27.CENPE", "TYMS", "IFI27_MX1", "IFI27_IFI44L")){
  tpm_sel <- get_zscore(variable=i, dataset=tpm_sel)
}


# For scores that go down in viral infection, multiply Z scores by x1----

tpm_sel$Herberg2 <- tpm_sel$Herberg2*-1
tpm_sel$Pennisi2 <- tpm_sel$Pennisi2*-1
tpm_sel$Lopez7 <- tpm_sel$Lopez7*-1
tpm_sel$Sampson10 <- tpm_sel$Sampson10*-1
tpm_sel$Li3 <- tpm_sel$Li3*-1
tpm_sel$Ravichandran10 <- tpm_sel$Ravichandran10*-1
tpm_sel$Rao8 <- tpm_sel$Rao8*-1
tpm_sel$Cappuccio11 <- tpm_sel$Cappuccio11*-1


# Re-create new signatures based on scaled data

tpm_sel <- tpm_sel %>%
  mutate(IFI27_MX1_scaled = (IFI27 + MX1)/2,
         IFI27_IFI44L_scaled = (IFI27 + IFI44L)/2)

for(i in c("IFI27_MX1_scaled", "IFI27_IFI44L_scaled")){
  tpm_sel <- get_zscore(variable=i, dataset=tpm_sel)
}


# Save signatures

saveRDS(tpm_sel, "tpm_sel.rds")
tpm_sel <- readRDS("tpm_sel.rds")

