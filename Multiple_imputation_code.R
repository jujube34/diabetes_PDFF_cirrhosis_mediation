################################################################ #
# Multiple imputation -----------------
################################################################ #
library(mice)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(magrittr)
library(regmedint)
library(broom)
library(glue)

# Variables
# DBbase: Presence of baseline diabetes.
# liver_mri_PDFF: Liver proton density fat fraction.
# f_Cirrhosis_sg: Follow-up time for incident cirrhosis.
# Cirrhosis_sgcase: Whether or not one developed cirrhosis in the follow-up.
# socioecoQ: socioeconomic status.
# Income: Annual household income.
# IPAQ: Physical activity.
# alc_g: Alcohol intake.
# smkcat: Smoking.
# HTNbase: Presence of baseline hypertension.
# Other variable names were explained in the Supplement of the article

# Read seed number and sensitivity parameter k
# from the command line interface
args <- commandArgs(T)
seed <- args[1] %>% as.numeric()
k <- args[2] %>% as.numeric()
if(is.na(k)) k <- 1
print("seed number: ")
print(seed)
set.seed(seed)

# Read data
load("20221224.RData")

array_outcome <- c("f_Cirrhosis_sg", "Cirrhosis_sgcase")

################################ #
# define function
################################ #
f.mediation.regmedint <- function(dt_mice_merged,
                                  avar = "DBbase",
                                  mvar = "liver_mri_PDFF",
                                  mediator_level_for_cde = 5,
                                  cvar, cvar_condition,
                                  use_original_scale_PDFF = TRUE,
                                  interaction = TRUE,
                                  mediator_model, outcome_model) {
    dt_mice_merged %>%
        ## Stacked up dataset
        mice::complete("long") %>%
        as_tibble() -> dt
    if (use_original_scale_PDFF == TRUE) {
        dt %<>% mutate(liver_mri_PDFF = exp(liver_mri_PDFF) - 1)
    }
    dt %>%
        group_by(.imp) %>%
        ## Nested data frame
        nest() %>%
        mutate(fit = map(data, function(data) {
            regmedint(
                data = data,
                ## Variables
                avar = avar,
                mvar = mvar,
                yvar = array_outcome[1],
                cvar = cvar,
                eventvar = array_outcome[2],

                ## Values at which effects are evaluated
                a0 = 0,
                a1 = 1,
                m_cde = mediator_level_for_cde,
                c_cond = cvar_condition,

                ## Model types
                mreg = mediator_model,
                yreg = outcome_model,

                ## Additional specification
                interaction = interaction,
                casecontrol = FALSE
            )
        })) %>%
        ## Extract point estimates and variance estimates
        mutate(
            coef_fit = map(fit, coef),
            vcov_fit = map(fit, vcov)
        ) -> result
    return(result)
}

################################ #


var_auxiliary_PDFF <- c(
 "pace_maker", "hearing_aid", "Cochlear_implant",
 "believed_safe_abdominal_MRI", "feeling_heart_racing",
 "feeling_tired", "asbestosis", "death_record_origin",
 "acceleration_175", "cancer_pain", "conc_VLDL",
 "conc_largeVLDL", "conc_largeLDL", "FC_in_VLDL",
 "conc_verylargeHDL"
)

var_factor <- c(
 "pace_maker", "hearing_aid", "Cochlear_implant",
 "believed_safe_abdominal_MRI", "feeling_heart_racing",
 "feeling_tired", "asbestosis", "cancer_pain"
)

var_impute <- c(
 "ALT", "AST", "GGT", "TB", "Urate", "WBC", "PLT",
 "BMI", "HbA1c", "TG", "CRP", "GRS_DM", "GRS_Cirrhosis_num",
 "age", "sex", "alc_g", "socioecoQ", "edu", "Income", 
 "smkcat", "IPAQ", "liver_mri_PDFF", "DBbase", "cancerbase",
 "HTNbase", "CHDbase",  "strokebase", "f_Cirrhosis_sg",
 "Cirrhosis_sgcase", var_auxiliary_PDFF
)

var_impute <- unique(var_impute)

dt_med <- dt.1 %>% dplyr::select(all_of(var_impute))

print("The dataset used for multiple imputation:")
print(dt_med)

################################ #
# construct predictor matrix -----------------
################################ #
pred <- quickpred(dt_med,
 mincor = 0.1,
 minpuc = 0,
 include = c("age","sex","socioecoQ","edu","Income"),
 exclude = "liver_mri_PDFF")

# use all auxiliary variabels for imputation of liver PDFF
pred["liver_mri_PDFF",var_auxiliary_PDFF] <- 1
pred["liver_mri_PDFF",var_auxiliary_PDFF]

print("predictorMatrix:")
print(pred)

################################ #
# Coding of categorical variables -----------------
################################ #

dt_med %<>%
  mutate(
    liver_mri_PDFF = liver_mri_PDFF %>% add(1) %>% log(),
    Income = case_when(
      Income == 9 ~ NA_character_,
      TRUE ~ Income %>% as.character()
    ) %>% as.integer(),
    edu = case_when(
      edu == 1 ~ "3",
      edu == 2 ~ "0",
      edu == 3 ~ "2",
      edu == 4 ~ "1",
      edu == 5 ~ "0",
      edu == 9 ~ NA_character_,
      TRUE ~ NA_character_,
    ) %>% as.integer(),
    IPAQ = case_when(
      IPAQ == 9 ~ NA_character_,
      TRUE ~ IPAQ %>% as.character()
    ) %>% as.integer(),
    socioecoQ = case_when(
      socioecoQ == 9 ~ NA_character_,
      TRUE ~ socioecoQ %>% as.character()
    ) %>% as.integer(),
    smkcat = case_when(
      smkcat == 9 ~ NA_character_,
      TRUE ~ smkcat %>% as.character()
    ) %>% as.integer(),
    GRS_DM = GRS_DM %>% scale() %>% as.numeric(),
    pace_maker = pace_maker %>% as.character() %>% as.integer(),
    hearing_aid = hearing_aid %>% as.character() %>% as.integer(),
    Cochlear_implant = Cochlear_implant %>% as.character() %>% as.integer(),
    HTNbase = HTNbase %>% as.character() %>% as.integer(),
    CHDbase = CHDbase %>% as.character() %>% as.integer(),
    believed_safe_abdominal_MRI = case_when(
      believed_safe_abdominal_MRI == 9 ~ NA_character_,
      TRUE ~ believed_safe_abdominal_MRI %>% as.character()
    ) %>% as.integer(),
    feeling_heart_racing = case_when(
      feeling_heart_racing == -601 ~ 1,
      feeling_heart_racing == -602 ~ 1,
      feeling_heart_racing == -600 ~ 0,
      TRUE ~ NA_real_,
    ),
    feeling_tired = case_when(
      feeling_tired == -601 ~ 1,
      feeling_tired == -602 ~ 1,
      feeling_tired == -600 ~ 0,
      TRUE ~ NA_real_,
    ),
    asbestosis = asbestosis %>% as.character() %>% as.integer(),
    death_record_origin = case_when(
      death_record_origin == "E/W" ~ 0,
      death_record_origin == "SCOT" ~ 1,
      TRUE ~ NA_real_
    ),
    cancer_pain = cancer_pain %>% as.character() %>% as.integer()
  )

f.turn.neg.to.NA <- function(x) {
  idx_na <- x %>% is.na()
  x[idx_na] <- -1
  x <- case_when(
    x < 0 ~ NA_real_,
    TRUE ~ x
  )
}

# turn factors into integers
dt_med %<>% mutate_at(var_factor, f.turn.neg.to.NA)

################################ #
# start multiple imputation -----------------
################################ #
# We make 5 imputations per thread of R program,
# and invoked 10 threads in total,
# which produced 50 imputed data sets.
m <- 5 
print("number of imputation:")
print(m)
T <- 20
print("number of iteration:")
print(T)

imp_pm <- mice::mice(
 dt_med,
 m = m,
 maxit = 0,
 predictorMatrix = pred,
 seed = seed,
 print = FALSE)

# post-processing
# When k = 1.0, no change was made to the imputed data set.
post <- imp_pm$post
post["liver_mri_PDFF"] <- paste("imp[[j]][,i] <-", k, "* imp[[j]][,i]")
imp_pm$post <- post

pnde <- matrix(NA, nrow = T, ncol = m)
tnie <- matrix(NA, nrow = T, ncol = m)
te <- matrix(NA, nrow = T, ncol = m)
pm <- matrix(NA, nrow = T, ncol = m)

# confounders of mediation effect 
cvar <- c("age", "sex")
cvar_condition <- c(50, 0)
array_outcome <- c("f_Cirrhosis_sg", "Cirrhosis_sgcase")

for (i in 1:T) {
    cat("Iteration = ", i,"/", T,"seed = ",seed,"\n")

    imp_pm <- mice::mice.mids(imp_pm, maxit = 1, seed = seed, print = FALSE)
    # calculate mediation effect of liver PDFF after each iteration to check for convergence
    mice_regmedint_Cox <- f.mediation.regmedint(
       dt_mice_merged = imp_pm,
        cvar = cvar,
        mvar = "liver_mri_PDFF",
        mediator_level_for_cde = 5,
        cvar_condition = cvar_condition,
        mediator_model = "linear",
        outcome_model = "survCox",
        interaction = TRUE, 
    )
    pnde[i, ] <- mice_regmedint_Cox$coef_fit %>% sapply(.,function(x)x["pnde"])
    tnie[i, ] <- mice_regmedint_Cox$coef_fit %>% sapply(.,function(x)x["tnie"])
    te[i, ] <- mice_regmedint_Cox$coef_fit %>% sapply(.,function(x)x["te"])
    pm[i, ] <- mice_regmedint_Cox$coef_fit %>% sapply(.,function(x)x["pm"])
}

print("Multiple imputation completed.")

################################ #
# save data -----------------
################################ #
time <- Sys.time() %>% format("%Y%m%d_%H%M%S")
filename <- glue("~/MI_data_{time}_seed_{seed}_k_{k}.RDS")
saveRDS(
    list(
        cde = cde,
        pnde = pnde,
        tnie = tnie,
        te = te,
        pm = pm,
        dt_mice_merged = imp_pm,
        k = k
    ),
    filename
)

