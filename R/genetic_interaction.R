library(DESeq2)
library(tidyverse)
library(here)
library(BART)

dds = readRDS(here("data/mimic_cc_dds_lrt.rds"))

cc_df = read_csv("../shared_data/Combined_Eds1_Rgt1_Lys14_CCTargets.csv")

shrunk_res_list = readRDS(here("data/shrunken_res_lists.rds"))

source(here("R/produce_results_tables.R"))

gi_res = extractGeneticInteractionTables(dds)

bart_df = gi_res$interaction_pattern %>%
  left_join(cc_df, by = c('genes' = 'id')) %>%
  left_join(as_tibble(shrunk_res_list$plus_lys$EDS1, rownames = 'genes')) %>%
  select(genes,pattern_rgt1, log2FoldChange, eds1_plus_lys_pval, rgt1_plus_lys_pval)

ggplot(bart_df) + geom_point(aes(-log(eds1_plus_lys_pval), log2FoldChange, color = pattern_rgt1))
#
# bart_df = bart_df[complete.cases(bart_df),]
#
# bart_out = BART::mc.mbart(
#   x.train = as.data.frame(select(bart_df, -pattern_rgt1)),
#   y.train = as.numeric(factor(pull(bart_df, pattern_rgt1), labels = seq(1,length(unique(pull(bart_df, pattern_rgt1)))))),
#   mc.cores = 10
# )
#
# # write_rds(bart_out, here("data/bart_pattern_onto_cc_eds_rgt1.rds"))
#
# ## load the advanced lung cancer example
# data(lung)
#
# group <- -which(is.na(lung[ , 7])) ## remove missing row for ph.karno
# times <- lung[group, 2]   ##lung$time
# delta <- lung[group, 3]-1 ##lung$status: 1=censored, 2=dead
# ##delta: 0=censored, 1=dead
#
# ## this study reports time in days rather than months like other studies
# ## coarsening from days to months will reduce the computational burden
# times <- ceiling(times/30)
#
# summary(times)
# table(delta)
#
# x.train <- as.matrix(lung[group, c(4, 5, 7)]) ## matrix of observed covariates
#
# ## lung$age:        Age in years
# ## lung$sex:        Male=1 Female=2
# ## lung$ph.karno:   Karnofsky performance score (dead=0:normal=100:by=10)
# ##                  rated by physician
#
# dimnames(x.train)[[2]] <- c('age(yr)', 'M(1):F(2)', 'ph.karno(0:100:10)')
#
# summary(x.train[ , 1])
# table(x.train[ , 2])
# table(x.train[ , 3])
#
# x.test <- matrix(nrow=84, ncol=3) ## matrix of covariate scenarios
#
# dimnames(x.test)[[2]] <- dimnames(x.train)[[2]]
#
# i <- 1
#
# for(age in 5*(9:15)) for(sex in 1:2) for(ph.karno in 10*(5:10)) {
#   x.test[i, ] <- c(age, sex, ph.karno)
#   i <- i+1
# }
#
# ## this x.test is relatively small, but often you will want to
# ## predict for a large x.test matrix which may cause problems
# ## due to consumption of RAM so we can predict separately
#
# ## mcparallel/mccollect do not exist on windows
# if(.Platform$OS.type=='unix') {
#   ##test BART with token run to ensure installation works
#   set.seed(99)
#   post <- surv.bart(x.train=x.train, times=times, delta=delta, nskip=5, ndpost=5, keepevery=1)
#
#   pre <- surv.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)
#
#   pred <- predict(post, pre$tx.test)
#   ##pred. <- surv.pwbart(pre$tx.test, post$treedraws, post$binaryOffset)
# }
#
# ## Not run:
# ## run one long MCMC chain in one process
# set.seed(99)
# post <- surv.bart(x.train=x.train, times=times, delta=delta)
#
# ## run "mc.cores" number of shorter MCMC chains in parallel processes
# ## post <- mc.surv.bart(x.train=x.train, times=times, delta=delta,
# ##                      mc.cores=5, seed=99)
#
# pre <- surv.pre.bart(x.train=x.train, times=times, delta=delta, x.test=x.test)
#
# pred <- predict(post, pre$tx.test)
#
# ## let's look at some survival curves
# ## first, a younger group with a healthier KPS
# ## age 50 with KPS=90: males and females
# ## males: row 17, females: row 23
# x.test[c(17, 23), ]
#
# low.risk.males <- 60*post$K+1:post$K ## K=unique times including censoring
# low.risk.females <- 6*post$K+1:post$K
#
# plot(post$times, pred$surv.test.mean[low.risk.males], type='s', col='blue',
#      main='Age 50 with KPS=90', xlab='t', ylab='S(t)', ylim=c(0, 1))
# points(post$times, pred$surv.test.mean[low.risk.females], type='s', col='red')
