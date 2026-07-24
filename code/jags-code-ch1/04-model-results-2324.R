# --------------------------------------------------------------------------------------------
##
## 04-model-results-2324
##
## Task: Summarize basic model results for OSS/ENES JAGS Occupancy Model 2023-2024
##
## Start Date: 02/26/2025
## Jasmine Williamson
##

## goals
# starting from scratch, using bits of Jay Jones' code

## insights
# model runs with 2023-2024 data
# made a parameter posteriors plot but i cant make sense of it

## settings --------------------------------------------------------------------------------

    rm(list=ls())
    library(R2jags)
    library(grid)
    library(gridExtra)
    library(abind)
    library(tidyr)
    library(ggplot2)


## load data --------------------------------------------------------------------------------
    
    load("data/jags/oss_enes_data_packaged_jw.RData") # my input data
    
    load("data/jags/jags-occu-model-objects.RData") # my model output


## extract results --------------------------------------------------------------------------------
    
    # extract outputs
    oh1.sl <- out.oh1$BUGSoutput$sims.list
    eh1.sl <- out.eh1$BUGSoutput$sims.list
    
    # parameters of interest
    coefficients <- out.oh1$BUGSoutput$summary[rownames(out.oh1$BUGSoutput$summary) %in% 
                                 c("beta0", "betaTFCL", "betaTFNC", "betaOWPB", "betaOWODF", "betaOWBLM", "betaYr24", "betaHU",
                                   "betaBU", "betaHB", "betaBS", "betaDW", "SalvEffect", "BurnEffect",
                                   "mu.a0", "mu.a1"),]
    
    # convert log-odds to probabilities
    logit_to_prob <- function(log_odds) {
      exp(log_odds)/(1 + exp(log_odds))
    }
    
    # calc median probabilities
    coeff <- as.data.frame(coefficients)
    coeff$probabilities <- logit_to_prob(coeff$mean)
    
    # plot
    df_probs <- data.frame(
      Parameter = rownames(coeff),
      Probability = coeff$probabilities,
      Lower = logit_to_prob(coeff$`2.5%`),
      Upper = logit_to_prob(coeff$`97.5%`)
    )
    
    # Sort by median probability
    df_probs <- df_probs[order(df_probs$Probability),]
    
    # Create the plot
    ggplot(df_probs, aes(x = reorder(Parameter, Probability), y = Probability)) +
      geom_pointrange(aes(ymin = Lower, ymax = Upper), color = "blue") +
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray") +
      coord_flip() +
      scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
      labs(title = "Posterior Distribution of Occupancy Probabilities",
           subtitle = "With 95% Credible Intervals",
           x = "Parameter",
           y = "Effect") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        