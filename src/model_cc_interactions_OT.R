# load packages
library(ggplot2)
library(data.table)
library(gridExtra)
library(stringr)
library(dplyr)

library(emmeans)
#library(statmod)
library(glmmTMB)
# examine color brewer palette:
library(RColorBrewer)
# View the PuRd palette
display.brewer.pal(9, "PuRd")

# set working dir
wkdir <- "~/Desktop/StewartLab_projects/Pierre_G_spinalcord_sc"
setwd(wkdir)

# read in week 1
data_dir= "weights/week1"
# loop through files to get weights
dirs <- list.dirs(path = data_dir, 
                    full.names = TRUE, recursive = FALSE)
dirs
# create empty list
data_list1 <- list()
# loop
for (d in dirs) {
  file <- list.files(path= d, pattern = "\\weight_matrix.txt$",full.names = TRUE)
  data <- read.table(file, header=TRUE)
  filename <- basename(file)
  name <- str_extract(filename, "\\d+wk_+\\d")
  print(name)
  x <- colnames(data)
  rownames(data) <- x
  data_list1[[name]] <- data
}
data_list1

# read in sham
data_dir= "weights/sham"
# loop through files to get weights
dirs2 <- list.dirs(path = data_dir, 
                  full.names = TRUE, recursive = FALSE)
dirs2
# create empty list
data_list2 <- list()
# loop
for (d in dirs2) {
  file <- list.files(path= d, pattern = "\\weight_matrix.txt$",full.names = TRUE)
  data <- read.table(file, header=TRUE)
  filename <- basename(file)
  name <- str_extract(filename, "sham_+\\d+\\d")
  print(name)
  x <- colnames(data)
  rownames(data) <- x
  data_list2[[name]] <- data
}
names(data_list2)

# read in nonrecovered
data_dir= "weights/nonrecovered"
NRdirs <- list.dirs(path = data_dir, 
                  full.names = TRUE, recursive = FALSE)
NRdirs
data_listNR <- list()

for (d in NRdirs) {
  NRdirs2 <- list.dirs(path = d, 
                       full.names = TRUE, recursive = FALSE)
  #print(NRdirs2)
  for (d2 in NRdirs2) {
    file <- list.files(path= d2, pattern = "\\weight_matrix.txt$", full.names = TRUE)
    data <- read.table(file, header=TRUE)
    filename <- basename(file)
    name <- str_extract(filename, "\\d+wkNR_+\\d+\\d")
    print(name)
    x <- colnames(data)
    rownames(data) <- x
    data_listNR[[name]] <- data
  }
}
names(data_listNR)

# read in recovered
data_dir= "weights/recovered"
Rdirs <- list.dirs(path = data_dir, 
                    full.names = TRUE, recursive = FALSE)
Rdirs
data_listR <- list()

for (d in Rdirs) {
  Rdirs2 <- list.dirs(path = d, 
                       full.names = TRUE, recursive = FALSE)
  #print(NRdirs2)
  for (d2 in Rdirs2) {
    file <- list.files(path= d2, pattern = "\\weight_matrix.txt$", full.names = TRUE)
    data <- read.table(file, header=TRUE)
    filename <- basename(file)
    name <- str_extract(filename, "\\d+wkR_+\\d+\\d")
    if(is.na(name)){
      name <- str_extract(filename, "\\d+wkR_+\\d")
    }
    print(name)
    x <- colnames(data)
    rownames(data) <- x
    data_listR[[name]] <- data
  }
}
names(data_listR)

# combine early weeks with NR
data_listNR <- c(data_listNR, data_list1, data_list2)
names(data_listNR)
# combine early with R
data_listR <- c(data_listR, data_list1, data_list2)
names(data_listR)

## loop to get all targets of a sender
# create empty dataframe
sender <- "T_Cells_and_Neutrophils"
sender <- "Fibroblasts"
df_T <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("source", "target", "treatment","timepoint","rep","value")
colnames(df_T) <- x
# loop through NR datalist to get all NR values
for(i in 1:length(data_listNR)){
  data <- data_listNR[[i]]
  # subset data
  data_subset <- data[sender,]
  # get name (timepoint)
  name <- names(data_listNR)[i]
  n2 <- str_split_i(name,"_",1)
  if(n2 == "sham"){
    timepoint <- "0wk"
  } else {
    timepoint <- str_extract(name, "\\d+wk")  # Extracts "1wk", "2wk", etc.
  }
  rep <- str_split_i(name,"_",2)
  # convert data to long form
  long_df_sub <- melt(setDT(data_subset), 
                      variable.name = "target")
  # add source and timepoint data
  long_df_sub$source <- sender
  long_df_sub$treatment <- "nonrecovered"
  long_df_sub$timepoint <- timepoint
  long_df_sub$rep <- rep
  df_T <- rbind(df_T,long_df_sub)
}
# loop through R dataset get all R values
for(i in 1:length(data_listR)){
  data <- data_listR[[i]]
  # subset data
  data_subset <- data[sender,]
  # get name (timepoint)
  name <- names(data_listR)[i]
  n2 <- str_split_i(name,"_",1)
  if(n2 == "sham"){
    timepoint <- "0wk"
  } else {
    timepoint <- str_extract(name, "\\d+wk")  # Extracts "1wk", "2wk", etc.
  }
  rep <- str_split_i(name,"_",2)
  # convert data to long form
  long_df_sub <- melt(setDT(data_subset), 
                      variable.name = "target")
  # add source and timepoint data
  long_df_sub$source <- sender
  long_df_sub$treatment <- "recovered"
  long_df_sub$timepoint <- timepoint
  long_df_sub$rep <- rep
  df_T <- rbind(df_T,long_df_sub)
}
df_T$target <- as.character(df_T$target)
df_T

# write
write.table(df_T, file="Tcellsource_interactionstrengths_reps_OT.txt", sep="\t", 
            quote=F, row.names = F)
write.table(df_T, file="Fibroblastssource_interactionstrengths_reps_OT.txt", sep="\t", 
            quote=F, row.names = F)
# read back in table
df_T <- read.table(file = "Tcell_interactions/Tcellsource_interactionstrengths_reps_OT.txt",
                   sep="\t", header = TRUE)
df_T <- read.table(file = "Fibroblast_interactions/Fibroblastssource_interactionstrengths_reps_OT.txt", 
                   sep="\t", header =TRUE)
# loop through for plots and linear model
sender <- "T_Cells_and_Neutrophils"
sender <- "Fibroblasts"
targets <- unique(df_T$target)
outdir = "Tcell_interactions/"
outdir = "Fibroblast_interactions/"
#target1= as.character(targets[[1]])
#df_Tsub <- subset(df_T, target==target1)
for(k in 1:length(targets)){
  target1= as.character(targets[[k]])
  print(target1)
  df_Tsub <- subset(df_T, target==target1)
  my_title <- paste0("Celltype communication source: ", sender,
                     " target: ", target1)
  
  # Initialize model type tracker
  model_type <- "tweedie"
  model_converged <- FALSE
  
  # Check for all-zero groups
  zero_check <- df_Tsub %>%
    group_by(treatment, timepoint) %>%
    summarize(all_zero = all(value == 0), .groups = "drop") %>%
    filter(all_zero) %>%
    nrow()
  
  # Attempt Tweedie first (unless all-zero groups exist)
  # Tweedie model at 1.5 is a GLM that models zeros as Poisson and positive values as Gamma 
  if(zero_check == 0) {
    model <- tryCatch({
      glmmTMB(value ~ treatment * timepoint,
              family = tweedie(link = "log"),
              data = df_Tsub)
      
    }, error = function(e) NULL)
    
    if(!is.null(model)) model_converged <- (model$fit$convergence == 0)
  }
  
  # Fallback to zero-inflated Gamma with offset if:
  # 1) Tweedie failed, or 2) all-zero groups exist
  if(zero_check > 0 || !model_converged) {
    df_fixed <- df_Tsub %>%
      mutate(value_adj = ifelse(value == 0, 1e-10, value))
    
    model <- glmmTMB(value_adj ~ treatment * timepoint,
                     ziformula = ~1,
                     family = Gamma(link = "log"),
                     data = df_fixed)
    model_type <- "gamma_zi_offset"
    model_converged <- (model$fit$convergence == 0)
  }
  summary(model)
  model_type
  model_converged
  # Add model type to output for reference
  model_summary <- capture.output({
    cat("MODEL TYPE:", model_type, "\n\n", "MODEL Convergence: ", model_converged,
        "\n\n")
    summary(model)
  })
  
  # Save outputs with model type in filename
  writeLines(model_summary, 
             paste0(outdir, sender, "-", target1, "_", model_type, "_model_out.txt"))
  
  # Get estimated marginal means
  emm <- emmeans(model, ~ treatment * timepoint, type = "response")
  # Convert to dataframe for plotting
  emm_df <- as.data.frame(emm)
  # write
  write.csv(emm_df, file= paste0(outdir, sender,"-",target1,"_emm.csv"),
              row.names = F)
  # get contrast means
  # Get the interaction contrasts
  contrasts <- emmeans(model, ~ treatment * timepoint)
  pairwise <- contrast(contrasts, interaction = "pairwise", adjust = "none")
  pw_df <- as.data.frame(pairwise)
  # write df
  write.csv(pw_df, file= paste0(outdir, sender,"-",target1,"_interactions_con.csv"),
            row.names = F)
  # Test the treatment effect at each timepoint
  simple_effects <- emmeans(model, pairwise ~ treatment | timepoint)
  em_df <- as.data.frame(simple_effects$emmeans)
  # look at contrasts for recovered-nonrecovered at each timepoint
  con_df <- as.data.frame(simple_effects$contrasts)
  # write df
  write.csv(con_df, file= paste0(outdir, sender,"-",target1,"_contrasts_R-NR.csv"),
            row.names = F)
  # line plot
  # Plot the difference between treatments across timepoints
  # plot with CI and SE
  p1 <- emmip(model, treatment ~ timepoint, CIs = FALSE, type = "response") +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width = 0.1) + 
    theme_minimal() + labs(title = str_wrap(my_title, 40), x="Time point", 
         y = "Predicted interaction strength") + 
    scale_color_brewer(palette="Paired")
    #scale_color_manual(values = brewer.pal(9, "PuRd")[4:6])
  
  # make pdf
  pdf(file=paste0(outdir, "cc_interaction_strength_emmeans_", sender,
                  "-", target1,"_noCIs.pdf"),
      height=4, width=5)
  print(p1)
  dev.off()
}

## get all interactions for specific target
target <- "T_Cells_and_Neutrophils"
target <- "Fibroblasts"
# create empty dataframe
df_T <- data.frame(matrix(ncol = 6, nrow = 0))
x <- c("source", "target", "treatment","timepoint","rep","value")
colnames(df_T) <- x

# loop through NR datalist to get all NR values
for(i in 1:length(data_listNR)){
  data <- data_listNR[[i]]
  name <- names(data_listNR)[i]
  # subset data
  data_subset <- as.data.frame(data[target])
  colnames(data_subset) <- "value"
  # get name (timepoint)
  n2 <- str_split_i(name,"_",1)
  if(n2 == "sham"){
    timepoint <- "0wk"
  } else {
    timepoint <- str_extract(name, "\\d+wk")  # Extracts "1wk", "2wk", etc.
  }
  rep <- str_split_i(name,"_",2)
  # convert rownames to source column
  data_subset <- tibble::rownames_to_column(data_subset, "source")
  # add target and timepoint data
  data_subset$target <- target
  data_subset$treatment <- "nonrecovered"
  data_subset$timepoint <- timepoint
  data_subset$rep <- rep
  df_T <- rbind(df_T,data_subset)
}
# loop through R dataset get all R values
for(i in 1:length(data_listR)){
  data <- data_listR[[i]]
  name <- names(data_listNR)[i]
  # subset data
  data_subset <- as.data.frame(data[target])
  colnames(data_subset) <- "value"
  # get name (timepoint)
  n2 <- str_split_i(name,"_",1)
  if(n2 == "sham"){
    timepoint <- "0wk"
  } else {
    timepoint <- str_extract(name, "\\d+wk")  # Extracts "1wk", "2wk", etc.
  }
  rep <- str_split_i(name,"_",2)
  # convert rownames to source column
  data_subset <- tibble::rownames_to_column(data_subset, "source")
  # add target and timepoint data
  data_subset$target <- target
  data_subset$treatment <- "recovered"
  data_subset$timepoint <- timepoint
  data_subset$rep <- rep
  df_T <- rbind(df_T,data_subset)
}
df_T$source <- as.character(df_T$source)
df_T
# write
write.table(df_T, file="Tcelltarget_interactionstrengths_OT_reps.txt", sep="\t", 
            quote=F, row.names = F)
write.table(df_T, file="Fibroblaststarget_interactionstrengths_OT_reps.txt", sep="\t", 
            quote=F, row.names = F)
# read back in
df_T <- read.table(file = "Tcell_interactions/Tcelltarget_interactionstrengths_OT_reps.txt",
                   sep="\t", header = TRUE)
df_T <- read.table(file = "Fibroblast_interactions/Fibroblaststarget_interactionstrengths_OT_reps.txt", 
                   sep="\t", header =TRUE)
# set target
target <- "T_Cells_and_Neutrophils"
target <- "Fibroblasts"
# loop through for plots
sources <- unique(df_T$source)
#source1 <- "Glycinergic_Neurons"
outdir = "Tcell_interactions/" 
outdir = "Fibroblast_interactions/"
for(k in 1:length(sources)){
  source1= as.character(sources[[k]])
  print(source1)
  #source1 <- "Microglial_Cells"
  df_Tsub <- subset(df_T, source==source1)
  my_title <- paste0("Celltype communication source: ", source1,
                     " target: ", target)
  
  # Initialize model type tracker
  model_type <- "tweedie"
  model_converged <- FALSE
  
  # Check for all-zero groups
  zero_check <- df_Tsub %>%
    group_by(treatment, timepoint) %>%
    summarize(all_zero = all(value == 0), .groups = "drop") %>%
  filter(all_zero) %>%
  nrow()

  # Attempt Tweedie first (unless all-zero groups exist)
  if(zero_check == 0) {
  model <- tryCatch({
    glmmTMB(value ~ treatment * timepoint,
            family = tweedie(link = "log"),
            data = df_Tsub)
    
    }, error = function(e) NULL)
  
    if(!is.null(model)) model_converged <- (model$fit$convergence == 0)
  }

  # Fallback to zero-inflated Gamma with offset if:
  # 1) Tweedie failed, or 2) all-zero groups exist
  if(zero_check > 0 || !model_converged) {
    df_fixed <- df_Tsub %>%
      mutate(value_adj = ifelse(value == 0, 1e-10, value))
  
    model <- glmmTMB(value_adj ~ treatment * timepoint,
                   ziformula = ~1,
                   family = Gamma(link = "log"),
                   data = df_fixed)
    model_type <- "gamma_zi_offset"
    model_converged <- (model$fit$convergence == 0)
  }
  summary(model)
  model_type
  model_converged
  # Add model type to output for reference
  model_summary <- capture.output({
    cat("MODEL TYPE:", model_type, "\n\n", "MODEL Convergence: ", model_converged,
        "\n\n")
    summary(model)
  })

  # Save outputs with model type in filename
  writeLines(model_summary, 
           paste0(outdir, source1, "-", target, "_", model_type, "_model_out.txt"))
  
  # Get estimated marginal means
  emm <- emmeans(model, ~ treatment * timepoint, type = "response")
  # Convert to dataframe for plotting
  emm_df <- as.data.frame(emm)
  # write
  write.csv(emm_df, file= paste0(outdir, source1,"-",target,"_emm.csv"),
            row.names = F)
  
  # get contrast means
  # Get the interaction contrasts
  contrasts <- emmeans(model, ~ treatment * timepoint)
  pairwise <- contrast(contrasts, interaction = "pairwise", adjust = "none")
  pw_df <- as.data.frame(pairwise)
  # write df
  write.csv(pw_df, file= paste0(outdir, source1,"-",target,"_interactions_con.csv"),
            row.names = F)
  
  # Test the treatment effect at each timepoint
  simple_effects <- emmeans(model, pairwise ~ treatment | timepoint)
  em_df <- as.data.frame(simple_effects$emmeans)
  # look at contrasts for recovered-nonrecovered at each timepoint
  con_df <- as.data.frame(simple_effects$contrasts)
  # write df
  write.csv(con_df, file= paste0(outdir, source1,"-",target,"_contrasts_R-NR.csv"),
            row.names = F)
  
  # line plot
  # Plot the difference between treatments across timepoints
  # plot with CI and SE
  p1 <- emmip(model, treatment ~ timepoint, CIs = FALSE, type = "response") +
    geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width = 0.1) + 
    theme_minimal() + labs(title = str_wrap(my_title, 40), x="Time point", 
                           y = "Predicted interaction strength") + 
    #scale_color_brewer(palette="Paired") +
    scale_color_manual(values = brewer.pal(9, "PuRd")[4:6])
  
  # make pdf
  pdf(file=paste0(outdir, "cc_interaction_strength_emmeans_", source1,
                  "-", target,"_noCI.pdf"),
      height=4, width=5)
  print(p1)
  dev.off()
}


################################ linear model ##################################
# tweedie
# var power of 1.5 because between poisson (1) and gamma (2) and this accepts zeros
model <- glmmTMB(value ~ treatment * timepoint, 
                 family = tweedie(link = "log"), 
                 data = df_Tsub)
# Fit a zero-inflated Gamma GLM
# this model fits a gamma model for the positive part and a logistic model for the probability of zero.
model <- glmmTMB(value ~ treatment*timepoint, ziformula = ~1, family = 
                   Gamma(link = "log"), data = df_Tsub)
# create offset for all zeros
df_fixed <- df_Tsub %>%
  mutate(
    value_adj = ifelse(value == 0, 1e-10, value)
  )

model <- glmmTMB(
  value_adj ~ treatment * timepoint,
  ziformula = ~1,  # Model excess zeros
  family = Gamma(link = "log"),
  data = df_fixed  # Using offset-adjusted data
)

# Summarize the model
summary(model)
output <- capture.output(summary(model), file=NULL,append=FALSE)
# write

# Get estimated marginal means
emm <- emmeans(model, ~ treatment * timepoint, type = "response")
# Convert to dataframe for plotting
emm_df <- as.data.frame(emm)
# write

# Plot with confidence intervals
ggplot(emm_df, aes(x = timepoint, y = response, color = treatment, group = treatment)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  labs(y = "Predicted Value", x = "Timepoint") +
  theme_minimal()

# get contrast means
# Get the interaction contrasts
contrasts <- emmeans(model, ~ treatment * timepoint)
pairwise <- contrast(contrasts, interaction = "pairwise", adjust = "none")
pw_df <- as.data.frame(pairwise)
# write df

# Test the treatment effect at each timepoint
simple_effects <- emmeans(model, pairwise ~ treatment | timepoint)
em_df <- as.data.frame(simple_effects$emmeans)
# look at contrasts for recovered-nonrecovered at each timepoint
con_df <- as.data.frame(simple_effects$contrasts)
# write df

# Plot the difference between treatments across timepoints
# plot with CI and SE
emmip(model, treatment ~ timepoint, CIs = TRUE, type = "response") +
  geom_errorbar(aes(ymin = yvar - SE, ymax = yvar + SE), width = 0.1) + theme_minimal()

# Calculate treatment effect (B - A) at each timepoint
#treat_effect <- contrast(contrasts, method = "pairwise", by = "timepoint")
#plot(treat_effect)

