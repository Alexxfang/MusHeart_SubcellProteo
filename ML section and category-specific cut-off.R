#####ML section and category-specific cut-off##########
# following scripts cover information for Supp Table 12-17

library(tidyverse)       # For data manipulation, string handling, and visualization
library(BiocParallel)    # For parallel processing in Bioconductor workflows
library(MSnbase)         # For mass spectrometry data handling
library(pRoloc)          # For subcellular localization analysis
library(STRINGdb)        # For protein-protein interaction analysis
library(pRolocExtra)     # For additional pRoloc functionality
library(subcellularvis)  # For GO-based subcellular localization annotations


# Set timeout option to handle long-running computations
options(timeout=3600)

# Detect available cores and register for parallel processing
a<- detectCores()
registerDoParallel(a-2)


# Load input datasets: protein expression, marker list, and metadata
my_data_import <- read.csv("norm_transform_lopit_18092023.csv")
my_marker <- read.csv("Maker_list_25112023_mpv2.csv")
my_metadata <- read.csv("metadata_29082023.csv")

# Combine the quantitative data with marker information
true_combine <- my_data_import %>%
  left_join(my_marker, by = c("Protein.Group" = "Protein.Group")) %>%
  distinct(Protein.Group, .keep_all = TRUE) %>%
  mutate(Markers = replace_na(Markers, "unknown"))%>%
  mutate(markers = Markers)


# Export relevant subsets of data for downstream use
write.csv(true_combine %>% dplyr::select(Protein.Group, contains(c("Aug", "Sep")),-contains(c("Sum."))), "f1.csv")
write.csv(true_combine %>% dplyr::select(Protein.Group, Markers, markers), "f2.csv")
write.csv(my_metadata %>% dplyr::filter(grepl("Aug|Sep", sampleNames)), "f3.csv")

# Organize data for the pRoloc pipeline (subcellular localization)
Control_trial <- readMSnSet(exprsFile = "f1.csv", featureDataFile = "f2.csv", phenoDataFile = "f3.csv", sep = ",")


#### SVM optimization and classification ####
paramssvm <- svmOptimisation(Control_trial, fcol = "Markers",
                          times = 100, xval = 5,
                          verbose = FALSE,seed = 409036654,
                          BPPARAM = BiocParallel::SnowParam(stop.on.error=TRUE))

# Option to reproduce MusHeart_SubcellProteo result
paramssvm <- readRDS("paramssvm.rds")

# Collect F1 score data from the SVM optimization results
F1<- paramssvm@f1Matrices
svm_F1_data_frame <- data.frame()
for (iter in 1:length(F1)) {
  aaa = as.data.frame(F1[[iter]])
  aaa$iteration = iter
  aaa$sigma = rownames(aaa)
  aab  = aaa%>%
    gather(key = "cost", value = "svm_F1_score", -c(sigma,iteration))
  svm_F1_data_frame = bind_rows(svm_F1_data_frame,aab)
}

# Summarize F1 scores for SVM by sigma and cost
svm_F1_data_frame_sum <- svm_F1_data_frame%>%
  group_by(sigma, cost)%>%
  summarise(mean = mean(svm_F1_score),
            sd = sd(svm_F1_score),
            med = median(svm_F1_score),
            max = max(svm_F1_score),
            min = min(svm_F1_score))

# Sigma == 0.01, cost ==4 are the optimal parameters
# Extract iterations and parameters pairs information
aab<- as.data.frame(getF1Scores(paramssvm))
aab$iter<- rownames(aab)
for_calculation <- aab%>%filter(sigma == 0.01 & cost ==4)%>%mutate(iter = as.numeric(iter))

# Extract confusion matrices and calculate categorical F1 scores using optimal parameters
CM<- paramssvm@cmMatrices
svm_F1_category <- data.frame()
for (iter in for_calculation$iter) {
  aaa <- as.data.frame(CM[[iter]])
  aaa$iteration <- iter
  bbb <- aaa%>%group_by(Prediction)%>%mutate(test_sum = sum(Freq))%>%ungroup(Prediction)%>%
    group_by(Reference)%>%mutate(actual_sum = sum(Freq))%>%ungroup(Reference)%>%
    filter(Prediction==Reference)%>%
    mutate(precision = Freq/test_sum,recall = Freq/actual_sum, F1 = 2*(precision*recall)/(precision+recall))
  svm_F1_category = bind_rows(svm_F1_category,bbb)
}

# Perform SVM classification with optimal parameters
svmres <- svmClassification(Control_trial, fcol = "Markers",
                            sigma = 0.01, cost = 4)

# Option to reproduce MusHeart_SubcellProteo result
svmres <- readRDS("svmres.rds")

# Get SVM classification dataframe
p1 <- getPredictions(svmres, fcol = "svm", mcol = "Markers")
svm_p1 <- fData(p1)


#### RF optimization and classification ####
paramsrf <- rfOptimisation(Control_trial, fcol = "Markers",
                          times = 100, xval = 5,
                          verbose = FALSE,seed = 409036654,
                          BPPARAM = BiocParallel::SnowParam(stop.on.error=TRUE))

# Option to reproduce MusHeart_SubcellProteo result
paramsrf <- readRDS("paramsrf.rds")

# Collect F1 score data from the RF optimization results
F1<- paramsrf@f1Matrices
aaa<- F1[[1]]
rf_F1_data_frame <- data.frame()
for (iter in 1:length(F1)) {
  aaa = as.data.frame(F1[[iter]])
  aaa$iteration = iter
  aab  = aaa%>%
    gather(key = "mtry", value = "rf_F1_score", -c(iteration))
  rf_F1_data_frame = bind_rows(rf_F1_data_frame,aab)
}

# Summarize F1 scores for RF by mtry
rf_F1_data_frame_sum <- rf_F1_data_frame%>%
  group_by(mtry)%>%
  summarise(mean = mean(rf_F1_score),
            sd = sd(rf_F1_score),
            med = median(rf_F1_score),
            max = max(rf_F1_score),
            min = min(rf_F1_score))

# mtry == 4 is optimal parameter
# Extract iterations and parameters pairs information
aab<- as.data.frame(getF1Scores(paramsrf))
aab$iter<- rownames(aab)
for_calculation <- aab%>%filter(mtry == 4)%>%mutate(iter = as.numeric(iter))

# Extract confusion matrices and calculate categorical F1 scores using optimal parameters
CM<- paramsrf@cmMatrices
rf_F1_category <- data.frame()
for (iter in for_calculation$iter) {
  aaa <- as.data.frame(CM[[iter]])
  aaa$iteration <- iter
  bbb <- aaa%>%group_by(Prediction)%>%mutate(test_sum = sum(Freq))%>%ungroup(Prediction)%>%
    group_by(Reference)%>%mutate(actual_sum = sum(Freq))%>%ungroup(Reference)%>%
    filter(Prediction==Reference)%>%
    mutate(precision = Freq/test_sum,recall = Freq/actual_sum, F1 = 2*(precision*recall)/(precision+recall))
  rf_F1_category = bind_rows(rf_F1_category,bbb)
}

# Perform RF classification with optimal parameters
rfres <- rfClassification(Control_trial, fcol = "Markers",
                          mtry = 4)

# Option to reproduce MusHeart_SubcellProteo result
rfres <- readRDS("rfres.rds")

# Get RF classification dataframe
p1 <- getPredictions(rfres, fcol = "rf", mcol = "Markers")
rf_p1 <- fData(p1)


#### XGB optimization and classification ####
paramsxgb <- xgboostOptimisation(Control_trial, fcol = "markers",max_depth = 4:7,
                                times = 100, xval = 5,nrounds =100,seed = 409036654,
                                verbose = FALSE)

# Option to reproduce MusHeart_SubcellProteo result
paramsxgb <- readRDS("paramsxgb.rds")

# Collect F1 score data from the XGB optimization results
F1<- paramsxgb@f1Matrices
aaa<-F1[[1]]
xgb_F1_data_frame <- data.frame()
for (iter in 1:length(F1)) {
  aaa = as.data.frame(F1[[iter]])
  aaa$iteration = iter
  aaa$gamma = rownames(aaa)
  aab  = aaa%>%
    gather(key = "max_depth", value = "xgb_F1_score", -c(gamma,iteration))
  xgb_F1_data_frame = bind_rows(xgb_F1_data_frame,aab)
}

# Summarize F1 scores for RF by gamma and max_depth
xgb_F1_data_frame_sum <- xgb_F1_data_frame%>%
  group_by(gamma, max_depth)%>%
  summarise(mean = mean(xgb_F1_score),
            sd = sd(xgb_F1_score),
            med = median(xgb_F1_score),
            max = max(xgb_F1_score),
            min = min(xgb_F1_score))

# gamma == 3.2, max_depth ==6 are the optimal parameters
# Extract iterations and parameters pairs information
aab<- as.data.frame(getF1Scores(paramsxgb))
aab$iter<- rownames(aab)
for_calculation <- aab%>%filter(gamma == 3.2 & max_depth ==6)%>%mutate(iter = as.numeric(iter))

# Extract confusion matrices and calculate categorical F1 scores using optimal parameters
CM<- paramsxgb@cmMatrices
xgb_F1_category <- data.frame()
for (iter in for_calculation$iter) {
  aaa <- as.data.frame(CM[[iter]])
  aaa$iteration <- iter
  bbb <- aaa%>%group_by(Prediction)%>%mutate(test_sum = sum(Freq))%>%ungroup(Prediction)%>%
    group_by(Reference)%>%mutate(actual_sum = sum(Freq))%>%ungroup(Reference)%>%
    filter(Prediction==Reference)%>%
    mutate(precision = Freq/test_sum,recall = Freq/actual_sum, F1 = 2*(precision*recall)/(precision+recall))
  xgb_F1_category = bind_rows(xgb_F1_category,bbb)
}

# Perform XGB classification with optimal parameters
xgbres <- xgboostClassification(Control_trial, fcol = "markers", 
                                gamma = 3.2, max_depth = 6, scores = "prediction")

# Option to reproduce MusHeart_SubcellProteo result
xgbres <- readRDS("xgbres.rds")

# Get XGB classification dataframe
p1 <- getPredictions(xgbres, fcol = "xgboost", mcol = "markers")
xgb_p1 <- fData(p1)

#### Combine all machine learning prediction results into a single dataframe ####
# Add Protein.Group to each model's results for merging
svm_p1$Protein.Group <-rownames(svm_p1)
rf_p1$Protein.Group <-rownames(rf_p1)
xgb_p1$Protein.Group <-rownames(xgb_p1)

# Convert results into data frames for compatibility with merge
svm_p1<- as.data.frame(svm_p1)%>%select(-markers, -Markers)
rf_p1 <- as.data.frame(rf_p1)%>%select(-markers, -Markers)
xgb_p1<- as.data.frame(xgb_p1)%>%select(-markers, -Markers)

# Merge the prediction results from SVM, RF, and XGBoost models by Protein.Group
ML_output_list <- list(svm_p1,rf_p1,xgb_p1)
ML_output_df<- ML_output_list%>%reduce(left_join, by='Protein.Group')
my_csv_protein_combine <- true_combine%>%left_join(ML_output_df, by = "Protein.Group")

#### ML and category-specific cut-off based on STRINGdb ####
# Initialize STRINGdb object for protein-protein interaction data
string_db <- STRINGdb$new(species=10090,score_threshold=0)

# Prepare the input data by separating protein identifiers (for STRINGdb mapping)
true_combine_uniprot_ID <- true_combine%>%mutate(Protein_Index = Protein.Group)%>%
  separate_rows(Protein_Index, sep = ";")%>%dplyr::select(Protein_Index)%>%
  distinct()%>%as.data.frame()

# Map proteins to STRINGdb identifiers
my_string_mapping <- string_db$map(true_combine_uniprot_ID, "Protein_Index",
                                   takeFirst=FALSE, removeUnmappedRows=TRUE,
                                   quiet=FALSE)
ID_index <- true_combine_uniprot_ID%>%full_join(my_string_mapping, by = "Protein_Index")%>%filter(!is.na(STRING_id))

# Retrieve protein-protein interactions from STRINGdb
inter <- string_db$get_interactions(my_string_mapping$STRING_id)

# Option to reproduce MusHeart_SubcellProteo result
inter <- readRDS("inter.rds")

# Clean and process the interaction data to calculate combined interaction scores
protein_interaction_pairs <- inter%>%distinct()%>%
  full_join(ID_index, by = c("from" = "STRING_id"), relationship = "many-to-many")%>%
  full_join(ID_index, by = c("to" = "STRING_id"), relationship = "many-to-many")%>%
  dplyr::rename("from.Protein_Index" = "Protein_Index.x","to.Protein_Index" = "Protein_Index.y")%>%
  mutate(Index = paste(from.Protein_Index,to.Protein_Index, sep = "~"))%>%
  filter(!is.na(from)&!is.na(to))%>%
  group_by(Index)%>%summarise(combined_score = mean(combined_score))%>%
  separate(Index, into = c("from.Protein.Group","to.Protein.Group"), sep = "~")

# Combine protein interaction data with original dataset
my_csv_protein_string <- true_combine%>%dplyr::select(Protein.Group)%>%mutate(Protein_index = Protein.Group)%>%
  separate_rows(Protein_index, sep = ";")%>%distinct()%>%
  left_join(protein_interaction_pairs, by = c("Protein_index" = "from.Protein.Group"),relationship = "many-to-many")%>%
  filter(!is.na(to.Protein.Group))

# Generate a summary of protein interactors for each Protein Group
my_string_ID <- my_csv_protein_string%>%dplyr::select(c(Protein.Group,to.Protein.Group))%>%
  group_by(Protein.Group)%>%summarise(Interactors = paste(to.Protein.Group, collapse = "|"))

# Category-specific probability cutoff (SVM)
# Organizing binary column of interactors found/not found from the same category
my_csv_protein_combine_cate_svm <- my_csv_protein_combine%>%
  group_by(svm.pred)%>%
  mutate(Protein_svm_cate = paste(Protein.Group, collapse = ";"))%>%ungroup(svm.pred)%>%
  left_join(my_string_ID, by = "Protein.Group")%>%
  mutate(interactor_found = ifelse(str_detect(Protein_svm_cate, Interactors), 1,0))%>%
  filter(!is.na(interactor_found))

# Category-specific FDR cut-off based on commonly protein interactors
# Generate permutation tests to estimate False Discovery Rate (FDR)
random_testing = data.frame()
for (subcell_location in unique(my_csv_protein_combine_cate_svm$svm.pred)) {
  non_marker = my_csv_protein_combine_cate_svm%>%filter(!grepl(subcell_location,svm.pred))
  marker = my_csv_protein_combine_cate_svm%>%filter(grepl(subcell_location,svm.pred))
  protein_cate = marker%>%dplyr::select(svm.pred, Protein_svm_cate)%>%distinct()
  Test_premuation <- bind_rows(replicate(50, non_marker %>% sample_n(length(marker$Protein.Ids)), simplify=F), .id="Obs")%>%
    mutate(new_interact = protein_cate$Protein_svm_cate)%>%
    mutate(interactor_found2 = ifelse(str_detect(new_interact, Interactors), 1,0))%>%
    summarise(Random = 1-sum(interactor_found2)/length(interactor_found2), Location = protein_cate$svm.pred)
  random_testing = rbind(random_testing,Test_premuation)
}

random_testing1<- random_testing%>%mutate(Random_005=0.05*(Random/mean(Random)), 
                                          Random_01=0.1*(Random/mean(Random)))

# Compute FDR and determine cutoffs
my_csv_string_svm <- my_csv_protein_combine_cate_svm%>%
  left_join(random_testing1, by = c("svm.pred" = "Location"))%>%
  group_by(svm.pred)%>%
  arrange(desc(svm.scores), .by_group = TRUE)%>%
  mutate(rank_svm = order(svm.scores, decreasing = T))%>%
  mutate(svm_sum = cumsum(interactor_found))%>%
  mutate(svm_fdr = 1-svm_sum/rank_svm)%>%
  mutate(svm_cutoff_005 = ifelse(svm_fdr<=Random_005,"Y","N"))%>%#filter(svm_cutoff == "Y")%>%
  mutate(svm_cutoff_01 = ifelse(svm_fdr<=Random_01,"Y","N"))%>%#filter(svm_cutoff == "Y")%>%
  summarise(svm_cut_005 = min(svm.scores[svm_cutoff_005=="Y"]),
            svm_cut_01 = min(svm.scores[svm_cutoff_01=="Y"]))

# Category-specific probability cutoff (RF)
# Organizing binary column of interactors found/not found from the same category
my_csv_protein_combine_cate_rf <- my_csv_protein_combine%>%filter(rf.pred != "unknown")%>%
  group_by(rf.pred)%>%
  mutate(Protein_rf_cate = paste(Protein.Group, collapse = ";"))%>%ungroup(rf.pred)%>%
  left_join(my_string_ID, by = "Protein.Group")%>%
  mutate(interactor_found = ifelse(str_detect(Protein_rf_cate, Interactors), 1,0))%>%
  filter(!is.na(interactor_found))

# Category-specific FDR cut-off based on commonly protein interactors
# Generate permutation tests to estimate False Discovery Rate (FDR)
random_testing = data.frame()
for (subcell_location in unique(my_csv_protein_combine_cate_rf$rf.pred)) {
  non_marker = my_csv_protein_combine_cate_rf%>%filter(!grepl(subcell_location,rf.pred))
  marker = my_csv_protein_combine_cate_rf%>%filter(grepl(subcell_location,rf.pred))
  protein_cate = marker%>%dplyr::select(rf.pred, Protein_rf_cate)%>%distinct()
  Test_premuation <- bind_rows(replicate(50, non_marker %>% sample_n(length(marker$Protein.Ids)), simplify=F), .id="Obs")%>%
    mutate(new_interact = protein_cate$Protein_rf_cate)%>%
    mutate(interactor_found2 = ifelse(str_detect(new_interact, Interactors), 1,0))%>%
    summarise(Random = 1-sum(interactor_found2)/length(interactor_found2), Location = protein_cate$rf.pred)
  random_testing = rbind(random_testing,Test_premuation)
}

random_testing1<- random_testing%>%mutate(Random_005=0.05*(Random/mean(Random)), 
                                          Random_01=0.1*(Random/mean(Random)))

# Compute FDR and determine cutoffs for RF
my_csv_string_rf <- my_csv_protein_combine_cate_rf%>%
  left_join(random_testing1, by = c("rf.pred" = "Location"))%>%
  group_by(rf.pred)%>%
  arrange(desc(rf.scores), .by_group = TRUE)%>%
  mutate(rank_rf = order(rf.scores, decreasing = T))%>%
  mutate(rf_sum = cumsum(interactor_found))%>%
  mutate(rf_fdr = 1-rf_sum/rank_rf)%>%
  mutate(rf_cutoff_005 = ifelse(rf_fdr<=Random_005,"Y","N"))%>%#filter(rf_cutoff == "Y")%>%
  mutate(rf_cutoff_01 = ifelse(rf_fdr<=Random_01,"Y","N"))%>%#filter(rf_cutoff == "Y")%>%
  summarise(rf_cut_005 = min(rf.scores[rf_cutoff_005=="Y"]),
            rf_cut_01 = min(rf.scores[rf_cutoff_01=="Y"]))

# Category-specific probability cutoff (XGB)
# Organizing binary column of interactors found/not found from the same category
my_csv_protein_combine_cate_xgb <- my_csv_protein_combine%>%
  group_by(xgboost.pred)%>%
  mutate(Protein_xgb_cate = paste(Protein.Group, collapse = ";"))%>%ungroup(xgboost.pred)%>%
  left_join(my_string_ID, by = "Protein.Group")%>%
  mutate(interactor_found = ifelse(str_detect(Protein_xgb_cate, Interactors), 1,0))%>%
  filter(!is.na(interactor_found))

# Category-specific FDR cut-off based on commonly protein interactors
# Generate permutation tests to estimate False Discovery Rate (FDR)
random_testing = data.frame()
for (subcell_location in unique(my_csv_protein_combine_cate_xgb$xgboost.pred)) {
  non_marker = my_csv_protein_combine_cate_xgb%>%filter(!grepl(subcell_location,xgboost.pred))
  marker = my_csv_protein_combine_cate_xgb%>%filter(grepl(subcell_location,xgboost.pred))
  protein_cate = marker%>%dplyr::select(xgboost.pred, Protein_xgb_cate)%>%distinct()
  Test_premuation <- bind_rows(replicate(50, non_marker %>% sample_n(length(marker$Protein.Ids)), simplify=F), .id="Obs")%>%
    mutate(new_interact = protein_cate$Protein_xgb_cate)%>%
    mutate(interactor_found2 = ifelse(str_detect(new_interact, Interactors), 1,0))%>%
    summarise(Random = 1-sum(interactor_found2)/length(interactor_found2), Location = protein_cate$xgboost.pred)
  random_testing = rbind(random_testing,Test_premuation)
}
random_testing1<- random_testing%>%mutate(Random_005=0.05*(Random/mean(Random)), 
                                          Random_01=0.1*(Random/mean(Random)))

# Compute FDR and determine cutoffs for XGB
my_csv_string_xgb <- my_csv_protein_combine_cate_xgb%>%
  left_join(random_testing1, by = c("xgboost.pred" = "Location"))%>%
  group_by(xgboost.pred)%>%
  arrange(desc(xgboost.scores), .by_group = TRUE)%>%
  mutate(rank_xgb = order(xgboost.scores, decreasing = T))%>%
  mutate(xgb_sum = cumsum(interactor_found))%>%
  mutate(xgb_fdr = 1-xgb_sum/rank_xgb)%>%
  mutate(xgb_cutoff_005 = ifelse(xgb_fdr<=Random_005,"Y","N"))%>%#filter(xgb_cutoff == "Y")%>%
  mutate(xgb_cutoff_01 = ifelse(xgb_fdr<=Random_01,"Y","N"))%>%#filter(xgb_cutoff == "Y")%>%
  summarise(xgb_cut_005 = min(xgboost.scores[xgb_cutoff_005=="Y"]),
            xgb_cut_01 = min(xgboost.scores[xgb_cutoff_01=="Y"]))

# Option to reproduce MusHeart_SubcellProteo result
my_csv_string_svm <- readRDS("my_csv_string_svm.rds")
my_csv_string_rf <- readRDS("my_csv_string_rf.rds")
my_csv_string_xgb <- readRDS("my_csv_string_xgb.rds")


#### ML and category-specific cut-off based on GO-based protein subcellular annotations ####
# Launch the subcellular visualization app (if needed, for interactive visualization)
subcellularapp()

# Prepare the protein data for GO-based annotations
for_subcellR <- my_csv_protein_combine%>%dplyr::select(Protein.Group,Markers)%>%mutate(Protein_index = Protein.Group)%>%
  separate_rows(Protein_index, sep = ";")%>%distinct()

# Option to reproduce MusHeart_SubcellProteo result
subcell_GOCC <- read_csv("CellVis_fullEnrichment_2023_12_16.csv")
subcell_GOCC_xp <-subcell_GOCC%>%separate_rows(Genes, sep = ",")%>%separate_rows(group, sep = ", ")
subccell_gocc <- for_subcellR%>%left_join(subcell_GOCC_xp,by = c("Protein_index" = "Genes"),relationship = "many-to-many")%>%
  filter(!is.na(n))

# Organize category-specific subcellular annotations 
subccell_gocc_summary <- subccell_gocc%>%filter(Markers != "unknown")%>%
  mutate(Colannotation = ifelse(!grepl("PX|CYTO|ST|RCX",Markers), "non_Cyto","Cyto"))%>%
  filter(Colannotation!="non_Cyto"|group!="Cytoplasm")%>%
  mutate(Colannotation = ifelse(!grepl("PX|NU|RCX",Markers), "non_Nu","Nu"))%>%
  filter(Colannotation!="non_Nu"|group!="Nucleus")%>%
  group_by((Markers))%>%mutate(nObs = n())%>%ungroup(Markers)%>%
  group_by(Markers,group)%>%summarise(count_cate = n(),Obs = nObs)%>%group_by(Markers)%>%
  mutate(Ratio = count_cate/Obs)%>%
  arrange(desc(Ratio), .by_group = TRUE)%>%distinct()%>%
  filter(Ratio>=0.04)%>%
  summarise(Location = paste(group, collapse = "|"))

# Summarize subcellular localization for each Protein.Group
subccell_protein_summary <- subccell_gocc%>%group_by(Protein.Group)%>%
  summarise(Location = paste(group, collapse = "|"))

# Merge the subcellular localization data with the original dataset
my_csv_protein_combine_go <- my_csv_protein_combine%>%left_join(subccell_protein_summary, by = "Protein.Group")%>%
  filter(!is.na(Location))

# Category-specific FDR cut-off based on GOCC-annotations (SVM) 
# Add GOCC annotations to the SVM prediction results
my_csv_protein_combine_go_markers <- my_csv_protein_combine_go%>%filter(Markers != "unknown")%>%
  left_join(subccell_gocc_summary, by = c("svm.pred" = "Markers"))%>%
  mutate(gocc_found = ifelse(str_detect(Location.x, Location.y), 1,0))%>%
  group_by(svm.pred)%>%
  arrange(desc(gocc_found), .by_group = TRUE)%>%
  mutate(rank_svm = order(svm.scores, decreasing = T))%>%
  mutate(svm_sum = cumsum(gocc_found))%>%
  mutate(svm_fdr = svm_sum/rank_svm)%>%
  filter(svm_fdr<1)

# Calculate FDR-based cutoff for SVM model
my_csv_protein_combine_go_svm <- my_csv_protein_combine_go%>%
  left_join(subccell_gocc_summary, by = c("svm.pred" = "Markers"))%>%
  mutate(gocc_found = ifelse(str_detect(Location.x, Location.y), 1,0))%>%
  group_by(svm.pred)%>%
  arrange(desc(svm.scores), .by_group = TRUE)%>%
  mutate(rank_svm = order(svm.scores, decreasing = T))%>%
  mutate(svm_sum = cumsum(gocc_found))%>%
  mutate(svm_fdr = 1-svm_sum/rank_svm)%>%
  mutate(svm_cutoff_005 = ifelse(svm_fdr<=0.05,"Y","N"))%>%
  mutate(svm_cutoff_01 = ifelse(svm_fdr<=0.1,"Y","N"))%>%
  summarise(svm_cut_005 = min(svm.scores[svm_cutoff_005=="Y"]),
            svm_cut_01 = min(svm.scores[svm_cutoff_01=="Y"]))


# Category-specific FDR cut-off based on GOCC-annotations (RF) 
# Repeat the FDR cutoff process for Random Forest (RF) predictions
my_csv_protein_combine_go_rf <- my_csv_protein_combine_go%>%
  filter(rf.pred!="unknown")%>%
  left_join(subccell_gocc_summary, by = c("rf.pred" = "Markers"))%>%
  mutate(gocc_found = ifelse(str_detect(Location.x, Location.y), 1,0))%>%
  group_by(rf.pred)%>%
  arrange(desc(rf.scores), .by_group = TRUE)%>%
  mutate(rank_rf = order(rf.scores, decreasing = T))%>%
  mutate(rf_sum = cumsum(gocc_found))%>%
  mutate(rf_fdr = 1-rf_sum/rank_rf)%>%
  mutate(rf_cutoff_005 = ifelse(rf_fdr<=0.05,"Y","N"))%>%
  mutate(rf_cutoff_01 = ifelse(rf_fdr<=0.1,"Y","N"))%>%
  summarise(rf_cut_005 = min(rf.scores[rf_cutoff_005=="Y"]),
            rf_cut_01 = min(rf.scores[rf_cutoff_01=="Y"]))

# Category-specific FDR cut-off based on GOCC-annotations (XGB)
# Repeat the FDR cutoff process for XGB predictions
my_csv_protein_combine_go_xgb <- my_csv_protein_combine_go%>%
  left_join(subccell_gocc_summary, by = c("xgboost.pred" = "Markers"))%>%
  mutate(gocc_found = ifelse(str_detect(Location.x, Location.y), 1,0))%>%
  group_by(xgboost.pred)%>%
  arrange(desc(xgboost.scores), .by_group = TRUE)%>%
  mutate(rank_xgb = order(xgboost.scores, decreasing = T))%>%
  mutate(xgb_sum = cumsum(gocc_found))%>%
  mutate(xgb_fdr = 1-xgb_sum/rank_xgb)%>%
  mutate(xgb_cutoff_005 = ifelse(xgb_fdr<=0.05,"Y","N"))%>%
  mutate(xgb_cutoff_01 = ifelse(xgb_fdr<=0.1,"Y","N"))%>%
  summarise(xgb_cut_005 = min(xgboost.scores[xgb_cutoff_005=="Y"]),
            xgb_cut_01 = min(xgboost.scores[xgb_cutoff_01=="Y"]))


#### Combine STRING- and GO-based FDR cut-off score per category for each algorithm #####
# Combine the STRINGdb and GOCC FDR cutoff results for each model
my_svm_combine <- my_csv_protein_combine_go_svm%>%full_join(my_csv_string_svm, by = "svm.pred")%>%
  mutate(across(where(is.numeric), ~ replace_na(., 1)))%>%
  rowwise()%>%
  mutate(svm_cut_005 = min(svm_cut_005.x,svm_cut_005.y),
         svm_cut_01 = min(svm_cut_01.x,svm_cut_01.y))%>%
  select(-contains(c(".x",".y")))

my_rf_combine <- my_csv_protein_combine_go_rf%>%full_join(my_csv_string_rf, by = "rf.pred")%>%
  mutate(across(where(is.numeric), ~ replace_na(., 1)))%>%
  rowwise()%>%
  mutate(rf_cut_005 = min(rf_cut_005.x,rf_cut_005.y),
         rf_cut_01 = min(rf_cut_01.x,rf_cut_01.y))%>%
  select(-contains(c(".x",".y")))

my_xgb_combine <- my_csv_protein_combine_go_xgb%>%full_join(my_csv_string_xgb, by = "xgboost.pred")%>%
  mutate(across(where(is.numeric), ~ replace_na(., 1)))%>%
  rowwise()%>%
  mutate(xgb_cut_005 = min(xgb_cut_005.x,xgb_cut_005.y),
         xgb_cut_01 = min(xgb_cut_01.x,xgb_cut_01.y))%>%
  select(-contains(c(".x",".y")))


# Combine categorical-specific cut-off for SVM, RF and XGB-based prediction
Categorical_cutoff_combine <- my_svm_combine%>%left_join(my_rf_combine, by = c("svm.pred" = "rf.pred"))%>%
  left_join(my_xgb_combine, by = c("svm.pred" = "xgboost.pred"))%>%rename(Category = svm.pred)%>%
  mutate(across(where(is.numeric), ~ round(., 3)))




