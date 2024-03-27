################################################################################
#####
#####  Bayesian Hierarchical & Latent Variable Modeling
#####  Init code & functions
#####  Johann-Friedrich Salzmann
#####
################################################################################

# Just some basic stuff I always use to code faster

# rm(list = ls(envir = globalenv()),envir = globalenv())

if (!require("pacman")) install.packages("pacman")
library(pacman)

p_load(magrittr)
p_load(ggplot2)
p_load(tidyverse)

p_load(cmdstanr)
p_load(bayesplot)
p_load(posterior)
p_load(tidybayes)
p_load(deSolve)
p_load(GGally)

"%cin%" = function(x,y){str_detect(y,x)}
"%xin%" = function(x,y){x %cin% y | y %cin% x}
"%ni%" = Negate("%in%")
"%nci%" = Negate("%cin%")
"%nxi%" = Negate("%xin%")
"%.%" = function(x,y){paste(x,y,sep = "")}
"%?%" = function(x, y) list(x = x, y = y)
"%:%" = function(xy, z) if(xy$x) xy$y else z
cNA = as.character(NA)
as.numeric.factor = . %>% as.numeric(levels(.))[.]



#### importCmdstanPosterior =============================================================

# function arguments --------------------------------------------------------------------
# posterior_file: path to the output posterior CSV file
# parameters: parameters to import posterior for
# incl_warmup: whether to include warmup draws (requires warmup to be stored in output CSV)

# function ------------------------------------------------------------------------------
importCmdstanPosterior <- function(posterior_files, parameters, incl_warmup = FALSE) {
  # collect meta information
  meta <- vroom::vroom(file = posterior_files,
                       delim = ",",
                       n_max = 44,
                       show_col_types=F)
  warmup_saved <- ifelse(as.integer(stringr::str_extract(meta[[1]][9], "[[:digit:]]+")) == 1L, TRUE, FALSE)
  if (warmup_saved) {
    warmup_draws <- stringr::str_extract(meta[[1]][8], "[[:digit:]]+") %>%
      as.integer()
  }
  sampling_draws <- stringr::str_extract(meta[[1]][7], "[[:digit:]]+") %>%
    as.integer()
  # collect posterior with warmup draws
  if (incl_warmup == TRUE) {
    if (warmup_saved == FALSE) {
      stop("Warmup draws not available.")
    }
    posterior <- vroom::vroom(posterior_files,
                              delim = ",",
                              comment = "# ",
                              n_max = warmup_draws+sampling_draws,
                              col_select = starts_with(parameters),
                              show_col_types=F)
    gc()
  } else {
    # in the presence of warmup draws, collect posterior without warmup draws
    if (warmup_saved == TRUE) {
      posterior <- vroom::vroom(posterior_files,
                                delim = ",",
                                comment = "# ",
                                n_max = warmup_draws+sampling_draws,
                                col_select = starts_with(parameters),
                                show_col_types=F) %>%
        slice(warmup_draws+1:(warmup_draws+sampling_draws))
      gc()
      # absent warmup draws, collect posterior without warmup draws      
      } else {
        posterior <- vroom::vroom(posterior_files,
                                  delim = ",",
                                  comment = "# ",
                                  n_max = sampling_draws,
                                  col_select = starts_with(parameters),
                                  show_col_types=F)  
        gc()
      }
    }
  return(posterior)
}
