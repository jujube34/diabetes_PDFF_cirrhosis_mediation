################################################################ #
# Mediation analysis
################################################################ #
library(regmedint)
library(dplyr)
library(stringr)
library(tidyr)
library(magrittr)
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

f.mediation.regmedint <- function(dt_mice_merged,
                                  avar = "DBbase", mvar = "liver_mri_PDFF",
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
# covariate parameter setting -----------------
################################ #
cvar_crude <- c("age", "sex")
cvar_condition_crude <- c(50, 0)
cvar_adjusted <- c(cvar_crude, "socioecoQ", "Income", "IPAQ", "alc_g", "smkcat", "HTNbase")
cvar_condition_adjusted <- c(cvar_condition_crude, 0, 3, 1, 0, 0, 0)
array_outcome <- c("f_Cirrhosis_sg", "Cirrhosis_sgcase")

################################ #
# load imputed data sets for each k
################################ #
for (k in c("1.0","1.1", "1.2", "1.3", "1.4", "1.5")) {
    RDS_file_dir <- glue("R_import/PDFF_mediation_MI_datasets/Kdata/k{k}/")
    MI_files <- dir(RDS_file_dir)
    MI_files <- MI_files %>%
        str_ends(".RDS") %>%
        MI_files[.]
    dt_mice_merged <- NULL
    set.seed(1123)
    for (i in 2:length(MI_files)) {
        if (i == 2) {
            dt_mice1 <- readRDS(paste0(RDS_file_dir, MI_files[i]))
            dt_mice2 <- readRDS(paste0(RDS_file_dir, MI_files[i - 1]))
            dt_mice_merged <- mice::ibind(dt_mice1[[6]], dt_mice2[[6]])
        } else {
            paste0(RDS_file_dir, MI_files[i]) %>% print()
            dt_mice <- readRDS(paste0(RDS_file_dir, MI_files[i]))
            dt_mice_merged <- mice::ibind(dt_mice_merged, dt_mice[[6]])
        }
    }

    ############################### #
    # Cox regression for mediator model -----------------
    ################################ #
    # MRI-PDFF
    # Crude analysis
    mice_regmedint_Cox_PDFF_crude <- f.mediation.regmedint(
        dt_mice_merged = dt_mice_merged,
        cvar = cvar_crude,
        mvar = "liver_mri_PDFF",
        mediator_level_for_cde = 0,
        cvar_condition = cvar_condition_crude,
        mediator_model = "linear",
        outcome_model = "survCox",
        interaction = TRUE,
        use_original_scale_PDFF = TRUE
    )

    regmedint_mi <- mitools::MIcombine(
        results = mice_regmedint_Cox_PDFF_crude$coef_fit,
        variances = mice_regmedint_Cox_PDFF_crude$vcov_fit
    )

    dt_PDFF_crude <- f.summary.results(regmedint_mi)

    # MRI-PDFF
    # adjusted analysis
    mice_regmedint_Cox_PDFF_adj <- f.mediation.regmedint(
        dt_mice_merged = dt_mice_merged,
        cvar = cvar_adjusted,
        mvar = "liver_mri_PDFF",
        mediator_level_for_cde = 0,
        cvar_condition = cvar_condition_adjusted,
        mediator_model = "linear",
        outcome_model = "survCox",
        interaction = TRUE,
        use_original_scale_PDFF = TRUE
    )

    regmedint_mi <- mitools::MIcombine(
        results = mice_regmedint_Cox_PDFF_adj$coef_fit,
        variances = mice_regmedint_Cox_PDFF_adj$vcov_fit
    )

    dt_PDFF_adj <- f.summary.results(regmedint_mi)

    # BMI
    # Crude analysis
    mice_regmedint_Cox_BMI_crude <- f.mediation.regmedint(
        dt_mice_merged = dt_mice_merged,
        cvar = cvar_crude,
        mvar = "BMI",
        mediator_level_for_cde = 22,
        cvar_condition = cvar_condition_crude,
        mediator_model = "linear",
        outcome_model = "survCox",
        interaction = TRUE
    )

    regmedint_mi <- mitools::MIcombine(
        results = mice_regmedint_Cox_BMI_crude$coef_fit,
        variances = mice_regmedint_Cox_BMI_crude$vcov_fit
    )

    dt_BMI_crude <- f.summary.results(regmedint_mi)

    # BMI
    # adjusted analysis
    mice_regmedint_Cox_BMI_adj <- f.mediation.regmedint(
        dt_mice_merged = dt_mice_merged,
        cvar = cvar_adjusted,
        mvar = "BMI",
        mediator_level_for_cde = 22,
        cvar_condition = cvar_condition_adjusted,
        mediator_model = "linear",
        outcome_model = "survCox",
        interaction = TRUE
    )

    regmedint_mi <- mitools::MIcombine(
        results = mice_regmedint_Cox_BMI_adj$coef_fit,
        variances = mice_regmedint_Cox_BMI_adj$vcov_fit
    )

    dt_BMI_adj <- f.summary.results(regmedint_mi)

    ################################ #
    # save results -----------------
    ################################ #
    saveRDS(
     list( dt_PDFF_crude, dt_PDFF_adj, dt_BMI_crude, dt_BMI_adj),
     glue("Mediation_effects_k{k}.RDS")
    )
}