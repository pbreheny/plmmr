l <- readRDS('/mnt/cleft/debug/problem-data.rds')
fit <- plmm_fit(l$prep,
                l$std_X_details,
                l$eta_star,
                l$penalty_factor,
                l$fbm_flag,
                l$penalty)
boxplot(fit$std_scale_beta[,100])
which.min(fit$std_scale_beta[,100])
table(l$prep$std_X[,3879])
