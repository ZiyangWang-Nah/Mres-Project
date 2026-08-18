# Store AUC results for each fold, per species 
results <- data.frame(
  fold              = 1:5,
  auc_glm_Macro     =NA,
  auc_glm_MacroCp   =NA,
  auc_glm_MicroCp   =NA,
  auc_glm_Micro     =NA,
  auc_glm_Lc        =NA,
  auc_glm_Lc_int_Macro = NA,
  auc_glm_Lc_int_MacroCp = NA,
  auc_glm_Lc_int_MicroCp = NA,
  auc_glm_Lc_int_Micro = NA,
  auc_glm_Lc_add_Macro = NA,
  auc_glm_Lc_add_MacroCp = NA,
  auc_glm_Lc_add_MicroCp = NA,
  auc_glm_Lc_add_Micro = NA,
  auc_maxent_Macro =NA,
  auc_maxent_MacroCp = NA,
  auc_maxent_MicroCp = NA,
  auc_maxent_Micro = NA,
  auc_maxent_Lc    = NA,
  auc_maxent_Lc_add_Macro = NA,
  auc_maxent_Lc_add_MacroCp = NA,
  auc_maxent_Lc_add_MicroCp = NA,
  auc_maxent_Lc_add_Micro = NA,
  auc_maxent_Lc_int_Macro = NA,
  auc_maxent_Lc_int_MacroCp = NA,
  auc_maxent_Lc_int_MicroCp = NA,
  auc_maxent_Lc_int_Micro = NA,
  auc_rf_Macro   = NA,
  auc_rf_MacroCp = NA,
  auc_rf_MicroCp = NA,
  auc_rf_Micro = NA,
  auc_rf_Lc   = NA,
  auc_rf_Lc_int_Macro = NA,
  auc_rf_Lc_int_MacroCp = NA,
  auc_rf_Lc_int_MicroCp = NA,
  auc_rf_Lc_int_Micro = NA
)
folds_list = readRDS("D:/Paper9/code/foldslist_Pitangus_sulphuratus.rds")
points_climate= fread("D:/Paper9/result/Model Input/final/Pitangus_sulphuratus_FULL_Macro_Micro_CP_merged.csv")

# Fit the GLM, Maxent and RF with Macro, Micro and Land cover variables, with 35 predictor * algorithm sets
for (k in 1:5) {
  
  cat("\n====== Fold", k, "/5 ======\n")
  
  train_idx <- folds_list[[k]][[1]]
  test_idx  <- folds_list[[k]][[2]]
  
  train_data <- points_climate[train_idx, ]
  test_data  <- points_climate[test_idx,  ]
  
  cat("训练集：", nrow(train_data), "行 | 测试集：", nrow(test_data), "行\n")
  
  # ── 1. GLM Macroclimate Only ──────────────────────────────────────
  glm_Macro <- glm(
    pa ~ Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Macro <- predict(glm_Macro, newdata = test_data, type = "response")
  results$auc_glm_Macro[k] <- roc(test_data$pa, pred_glm_Macro, quiet = TRUE)$auc
  
  # ── 2. GLM MacroCp Only ──────────────────────────────────────
  glm_MacroCp <- glm(
    pa ~ Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmin + Macro_CP_SMmax,
    data   = train_data,
    family = binomial()
  )
  pred_glm_MacroCp <- predict(glm_MacroCp, newdata = test_data, type = "response")
  results$auc_glm_MacroCp[k] <- roc(test_data$pa, pred_glm_MacroCp, quiet = TRUE)$auc
  
  # ── 3. GLM MicroCp Only ──────────────────────────────────────
  glm_MicroCp <- glm(
    pa ~ Micro_CP_Tmin + Micro_CP_Tmax + Micro_CP_SMmin + Micro_CP_SMmax,
    data   = train_data,
    family = binomial()
  )
  pred_glm_MicroCp <- predict(glm_MicroCp, newdata = test_data, type = "response")
  results$auc_glm_MicroCp[k] <- roc(test_data$pa, pred_glm_MicroCp, quiet = TRUE)$auc
  
  # ── 4. GLM Microclimate Only ──────────────────────────────────────
  glm_Micro <- glm(
    pa ~ Micro_Bio5 + Micro_Bio6 + Micro_Bio13 + Micro_Bio14,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Micro <- predict(glm_Micro, newdata = test_data, type = "response")
  results$auc_glm_Micro[k] <- roc(test_data$pa, pred_glm_Micro, quiet = TRUE)$auc
  
  # ── 5. GLM Landcover Only ──────────────────────────────────────
  glm_Lc <- glm(
    pa ~ landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc <- predict(glm_Lc, newdata = test_data, type = "response")
  results$auc_glm_Lc[k] <- roc(test_data$pa, pred_glm_Lc, quiet = TRUE)$auc
  
  # ── 6. GLM Macroclimate * Landcover ──────────────────────────────────────
  glm_Lc_int_Macro <- glm(
    pa ~ (Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14) * landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_int_Macro  <- predict(glm_Lc_int_Macro, newdata = test_data, type = "response")
  results$auc_glm_Lc_int_Macro[k] <- roc(test_data$pa, pred_glm_Lc_int_Macro, quiet = TRUE)$auc
  
  # ── 7. GLM MacroCp * Landcover ──────────────────────────────────────
  glm_Lc_int_MacroCp <- glm(
    pa ~ (Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmin + Macro_CP_SMmax) * landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_int_MacroCp <- predict(glm_Lc_int_MacroCp, newdata = test_data, type = "response")
  results$auc_glm_Lc_int_MacroCp[k] <- roc(test_data$pa, pred_glm_Lc_int_MacroCp, quiet = TRUE)$auc
  
  # ── 8. GLM MicroCp * Landcover ──────────────────────────────────────
  glm_Lc_int_MicroCp <- glm(
    pa ~ (Micro_CP_Tmin + Micro_CP_Tmax + Micro_CP_SMmin + Micro_CP_SMmax) * landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_int_MicroCp <- predict(glm_Lc_int_MicroCp, newdata = test_data, type = "response")
  results$auc_glm_Lc_int_MicroCp[k] <- roc(test_data$pa, pred_glm_Lc_int_MicroCp, quiet = TRUE)$auc
  
  # ── 9. GLM Microclimate * Landcover ──────────────────────────────────────
  glm_Lc_int_Micro <- glm(
    pa ~ (Micro_Bio5 + Micro_Bio6 + Micro_Bio13 + Micro_Bio14) * landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_int_Micro  <- predict(glm_Lc_int_Micro, newdata = test_data, type = "response")
  results$auc_glm_Lc_int_Micro[k] <- roc(test_data$pa, pred_glm_Lc_int_Micro, quiet = TRUE)$auc
  
  # ── 10. GLM Macroclimate + Landcover ──────────────────────────────────────
  glm_Lc_add_Macro <- glm(
    pa ~ Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14 + landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_add_Macro  <- predict(glm_Lc_add_Macro, newdata = test_data, type = "response")
  results$auc_glm_Lc_add_Macro[k] <- roc(test_data$pa, pred_glm_Lc_add_Macro, quiet = TRUE)$auc
  
  # ── 11. GLM MacroCp + Landcover ──────────────────────────────────────
  glm_Lc_add_MacroCp <- glm(
    pa ~ Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmin + Macro_CP_SMmax + landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_add_MacroCp <- predict(glm_Lc_add_MacroCp, newdata = test_data, type = "response")
  results$auc_glm_Lc_add_MacroCp[k] <- roc(test_data$pa, pred_glm_Lc_add_MacroCp, quiet = TRUE)$auc
  
  # ── 12. GLM MicroCp + Landcover ──────────────────────────────────────
  glm_Lc_add_MicroCp <- glm(
    pa ~ Micro_CP_Tmin + Micro_CP_Tmax + Micro_CP_SMmin + Micro_CP_SMmax + landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_add_MicroCp <- predict(glm_Lc_add_MicroCp, newdata = test_data, type = "response")
  results$auc_glm_Lc_add_MicroCp[k] <- roc(test_data$pa, pred_glm_Lc_add_MicroCp, quiet = TRUE)$auc
  
  # ── 13. GLM Microclimate + Landcover ──────────────────────────────────────
  glm_Lc_add_Micro <- glm(
    pa ~ Micro_Bio5 + Micro_Bio6 + Micro_Bio13 + Micro_Bio14 + landcover,
    data   = train_data,
    family = binomial()
  )
  pred_glm_Lc_add_Micro  <- predict(glm_Lc_add_Micro, newdata = test_data, type = "response")
  results$auc_glm_Lc_add_Micro[k] <- roc(test_data$pa, pred_glm_Lc_add_Micro, quiet = TRUE)$auc
  
  # ── Maxent（用dismo，和之前一致）──────────────────────
  # Maxent需要单独的presence坐标和环境变量矩阵
  train_macro_maxent_p <- train_data[train_data$pa == 1, 
                               c("Macro_Bio5","Macro_Bio6","Macro_Bio13","Macro_Bio14","landcover")]
  train_macro_maxent_b  <- train_data[train_data$pa == 0, 
                               c("Macro_Bio5","Macro_Bio6","Macro_Bio13","Macro_Bio14","landcover")]
  test_macro_maxent       <- test_data[, c("Macro_Bio5","Macro_Bio6","Macro_Bio13","Macro_Bio14","landcover")]
  
  train_macrocp_maxent_p <- train_data[train_data$pa == 1, 
                                     c("Macro_CP_Tmin","Macro_CP_Tmax","Macro_CP_SMmin", "Macro_CP_SMmax","landcover")]
  train_macrocp_maxent_b  <- train_data[train_data$pa == 0, 
                                      c("Macro_CP_Tmin","Macro_CP_Tmax","Macro_CP_SMmin", "Macro_CP_SMmax","landcover")]
  test_macrocp_maxent       <- test_data[, c("Macro_CP_Tmin","Macro_CP_Tmax","Macro_CP_SMmin", "Macro_CP_SMmax","landcover")]
  
  train_microcp_maxent_p <- train_data[train_data$pa == 1, 
                                       c("Micro_CP_Tmin","Micro_CP_Tmax","Micro_CP_SMmin", "Micro_CP_SMmax","landcover")]
  train_microcp_maxent_b  <- train_data[train_data$pa == 0, 
                                        c("Micro_CP_Tmin","Micro_CP_Tmax","Micro_CP_SMmin", "Micro_CP_SMmax","landcover")]
  test_microcp_maxent       <- test_data[, c("Micro_CP_Tmin","Micro_CP_Tmax","Micro_CP_SMmin", "Micro_CP_SMmax","landcover")]
  
  train_micro_maxent_p <- train_data[train_data$pa == 1, 
                                     c("Micro_Bio5","Micro_Bio6","Micro_Bio13","Micro_Bio14","landcover")]
  train_micro_maxent_b  <- train_data[train_data$pa == 0, 
                                      c("Micro_Bio5","Micro_Bio6","Micro_Bio13","Micro_Bio14","landcover")]
  test_micro_maxent       <- test_data[, c("Micro_Bio5","Micro_Bio6","Micro_Bio13","Micro_Bio14","landcover")]
  
  # ── 14. maxent Macroclimate only ──────────────────────────────────────
  tryCatch({
    maxent_Macro <- maxent(
      x       = rbind(train_macro_maxent_p, train_macro_maxent_b)[, -5],  # 只用连续变量
      p       = c(rep(1, nrow(train_macro_maxent_p)), 
                  rep(0, nrow(train_macro_maxent_b))),
      factors = NULL,
      args    = c("randomtestpoints=0", "betamultiplier=1")
    )
    pred_maxent_Macro <- predict(maxent_Macro, test_macro_maxent[, -5])
    results$auc_maxent_Macro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Macro), 
                                    quiet = TRUE)$auc
 error = function(e) {
    cat("Maxent fold", k, "出错：", e$message, "\n")
  }})
  
  # ── 15. maxent MacroCp  ──────────────────────────────────────
  tryCatch({
    maxent_MacroCp <- maxent(
      x       = rbind(train_macrocp_maxent_p, train_macrocp_maxent_b)[, -5],  # 只用连续变量
      p       = c(rep(1, nrow(train_macrocp_maxent_p)), 
                  rep(0, nrow(train_macrocp_maxent_b))),
      factors = NULL,
      args    = c("randomtestpoints=0", "betamultiplier=1")
    )
    pred_maxent_MacroCp <- predict(maxent_MacroCp, test_macrocp_maxent[, -5])
    results$auc_maxent_MacroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_MacroCp), 
                                       quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 16. maxent MicroCp  ──────────────────────────────────────
  tryCatch({
    maxent_MicroCp <- maxent(
      x       = rbind(train_microcp_maxent_p, train_microcp_maxent_b)[, -5],  # 只用连续变量
      p       = c(rep(1, nrow(train_microcp_maxent_p)), 
                  rep(0, nrow(train_microcp_maxent_b))),
      factors = NULL,
      args    = c("randomtestpoints=0", "betamultiplier=1")
    )
    pred_maxent_MicroCp <- predict(maxent_MicroCp, test_microcp_maxent[, -5])
    results$auc_maxent_MicroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_MicroCp), 
                                         quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 17. maxent Microclimate only ──────────────────────────────────────
  tryCatch({
    maxent_Micro <- maxent(
      x       = rbind(train_micro_maxent_p, train_micro_maxent_b)[, -5],  # 只用连续变量
      p       = c(rep(1, nrow(train_micro_maxent_p)), 
                  rep(0, nrow(train_micro_maxent_b))),
      factors = NULL,
      args    = c("randomtestpoints=0", "betamultiplier=1")
    )
    pred_maxent_Micro <- predict(maxent_Micro, test_micro_maxent[, -5])
    results$auc_maxent_Micro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Micro), 
                                       quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 18. landcover only ──────────────────────────────────────
  tryCatch({
    lc_train <- data.frame(
      landcover = as.factor(c(train_micro_maxent_p$landcover, 
                              train_micro_maxent_b$landcover))
    )
    lc_test  <- data.frame(
      landcover = as.factor(test_micro_maxent$landcover)
    )
    
    maxent_lc <- maxent(
      x       = lc_train,
      p       = c(rep(1, nrow(train_micro_maxent_p)), 
                  rep(0, nrow(train_micro_maxent_b))),
      factors = "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1")
    )
    pred_maxent_lc <- predict(maxent_lc, lc_test)
    results$auc_maxent_Lc[k] <- roc(test_data$pa, as.numeric(pred_maxent_lc), 
                                    quiet = TRUE)$auc
    error = function(e) {
    cat("Maxent fold", k, "出错：", e$message, "\n")
  }})
  
  # ── 19. maxent Macroclimate + Landcover ──────────────────────────────────────
  tryCatch({
    maxent_Lc_add_Macro <- maxent(
      x       = rbind(train_macro_maxent_p, train_macro_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_macro_maxent_p)), 
                  rep(0, nrow(train_macro_maxent_b))),
      factors =  "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=false",        # 禁用 Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_add_Macro <- predict(maxent_Lc_add_Macro, test_macro_maxent)
    results$auc_maxent_Lc_add_Macro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_add_Macro), 
                                       quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 20. maxent MacroCp + Landcover  ──────────────────────────────────────
  tryCatch({
    maxent_Lc_add_MacroCp <- maxent(
      x       = rbind(train_macrocp_maxent_p, train_macrocp_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_macrocp_maxent_p)), 
                  rep(0, nrow(train_macrocp_maxent_b))),
      factors = "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=false",        # 禁用 Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_add_MacroCp <- predict(maxent_Lc_add_MacroCp, test_macrocp_maxent)
    results$auc_maxent_Lc_add_MacroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_add_MacroCp), 
                                         quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 21. maxent MicroCp + landcover  ──────────────────────────────────────
  tryCatch({
    maxent_Lc_add_MicroCp <- maxent(
      x       = rbind(train_microcp_maxent_p, train_microcp_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_microcp_maxent_p)), 
                  rep(0, nrow(train_microcp_maxent_b))),
      factors = "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=false",        # 禁用 Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_add_MicroCp <- predict(maxent_Lc_add_MicroCp, test_microcp_maxent)
    results$auc_maxent_Lc_add_MicroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_add_MicroCp), 
                                         quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 22. maxent Microclimate + landcover ──────────────────────────────────────
  tryCatch({
    maxent_Lc_add_Micro <- maxent(
      x       = rbind(train_micro_maxent_p, train_micro_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_micro_maxent_p)), 
                  rep(0, nrow(train_micro_maxent_b))),
      factors =  "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=false",        # 禁用 Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_add_Micro <- predict(maxent_Lc_add_Micro, test_micro_maxent)
    results$auc_maxent_Lc_add_Micro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_add_Micro), 
                                       quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 23. maxent Macroclimate * Landcover ──────────────────────────────────────
  tryCatch({
    maxent_Lc_int_Macro <- maxent(
      x       = rbind(train_macro_maxent_p, train_macro_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_macro_maxent_p)), 
                  rep(0, nrow(train_macro_maxent_b))),
      factors =  "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=true",        # turn on Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_int_Macro <- predict(maxent_Lc_int_Macro, test_macro_maxent)
    results$auc_maxent_Lc_int_Macro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_int_Macro), 
                                              quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 24. maxent MacroCp * Landcover  ──────────────────────────────────────
  tryCatch({
    maxent_Lc_int_MacroCp <- maxent(
      x       = rbind(train_macrocp_maxent_p, train_macrocp_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_macrocp_maxent_p)), 
                  rep(0, nrow(train_macrocp_maxent_b))),
      factors = "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=true",        # turn on Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_int_MacroCp <- predict(maxent_Lc_int_MacroCp, test_macrocp_maxent)
    results$auc_maxent_Lc_int_MacroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_int_MacroCp), 
                                                quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 25. maxent MicroCp * landcover  ──────────────────────────────────────
  tryCatch({
    maxent_Lc_int_MicroCp <- maxent(
      x       = rbind(train_microcp_maxent_p, train_microcp_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_microcp_maxent_p)), 
                  rep(0, nrow(train_microcp_maxent_b))),
      factors = "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=true",        # turn on Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_int_MicroCp <- predict(maxent_Lc_int_MicroCp, test_microcp_maxent)
    results$auc_maxent_Lc_int_MicroCp[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_int_MicroCp), 
                                                quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── 26. maxent Microclimate * landcover ──────────────────────────────────────
  tryCatch({
    maxent_Lc_int_Micro <- maxent(
      x       = rbind(train_micro_maxent_p, train_micro_maxent_b),  # 只用连续变量
      p       = c(rep(1, nrow(train_micro_maxent_p)), 
                  rep(0, nrow(train_micro_maxent_b))),
      factors =  "landcover",
      args    = c("randomtestpoints=0", "betamultiplier=1",
                  "autofeature=false",    # 关闭自动选择，改成手动控制
                  "linear=true",          # 开启线性特征
                  "quadratic=true",       # 开启二次方特征
                  "hinge=true",           # 开启铰链特征
                  "product=true",        # turn on Product (Cross-product/交互) 特征！
                  "threshold=false")
    )
    pred_maxent_Lc_int_Micro <- predict(maxent_Lc_int_Micro, test_micro_maxent)
    results$auc_maxent_Lc_int_Micro[k] <- roc(test_data$pa, as.numeric(pred_maxent_Lc_int_Micro), 
                                              quiet = TRUE)$auc
    error = function(e) {
      cat("Maxent fold", k, "出错：", e$message, "\n")
    }})
  
  # ── Random Forest（下采样版本）────────────────────────
  n_balance <- min(sum(train_data$pa == 1), 
                   sum(train_data$pa == 0))  # 取较小的那个
  
  bg_idx_down   <- sample(which(train_data$pa == 0), size = n_balance, replace = FALSE)
  pres_idx_down <- sample(which(train_data$pa == 1), size = n_balance, replace = FALSE)
  
  train_down <- train_data[c(pres_idx_down, bg_idx_down), ]
  train_down$pa_factor <- as.factor(
    ifelse(train_down$pa == 1, "presence", "background")
  )
  
  cat("下采样后训练集：", nrow(train_down),
      "行（presence:", n_balance, "，bg:", n_balance, "）\n")
  
  # ── 27. RF Macroclimate only ──────────────────────────────────────
  rf_Macro <- randomForest(
    pa_factor ~ Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14,
    data       = train_down,
    ntree      = 500,
    mtry       = 2,         
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Macro <- predict(rf_Macro, test_data, type = "prob")[, "presence"]
  results$auc_rf_Macro[k] <- roc(test_data$pa, pred_rf_Macro, 
                              quiet = TRUE)$auc
  
  # ── 28. RF MacroCp only ──────────────────────────────────────
  rf_MacroCp <- randomForest(
    pa_factor ~ Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmax + Macro_CP_SMmin,
    data       = train_down,
    ntree      = 500,
    mtry       = 2,         
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_MacroCp <- predict(rf_MacroCp, test_data, type = "prob")[, "presence"]
  results$auc_rf_MacroCp[k] <- roc(test_data$pa, pred_rf_MacroCp, 
                                 quiet = TRUE)$auc
  
  # ── 29. RF MicroCp only ──────────────────────────────────────
  rf_MicroCp <- randomForest(
    pa_factor ~ Micro_CP_Tmin + Micro_CP_Tmax + Micro_CP_SMmin + Micro_CP_SMmax,
    data       = train_down,
    ntree      = 500,
    mtry       = 2,         
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_MicroCp <- predict(rf_MicroCp, test_data, type = "prob")[, "presence"]
  results$auc_rf_MicroCp[k] <- roc(test_data$pa, pred_rf_MicroCp, 
                                   quiet = TRUE)$auc
  # ── 30. RF Microclimate only ──────────────────────────────────────
  rf_Micro <- randomForest(
    pa_factor ~ Micro_Bio5 + Micro_Bio6 + Micro_Bio13 + Micro_Bio14,
    data       = train_down,
    ntree      = 500,
    mtry       = 2,         
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Micro <- predict(rf_Micro, test_data, type = "prob")[, "presence"]
  results$auc_rf_Micro[k] <- roc(test_data$pa, pred_rf_Micro, 
                                 quiet = TRUE)$auc
  
  # ── 31. RF Landcover only ──────────────────────────────────
  rf_Lc <- randomForest(
    pa_factor ~ landcover,
    data       = train_down,
    ntree      = 500,
    mtry       = 1,          # 只有1个变量
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Lc <- predict(rf_Lc, test_data, type = "prob")[, "presence"]
  results$auc_rf_Lc[k] <- roc(test_data$pa, pred_rf_Lc, 
                              quiet = TRUE)$auc
  
  # ── 32. RF Macroclimate * landcover ──────────────────────────────────────
  rf_Lc_int_Macro <- randomForest(
    pa_factor ~ Macro_Bio5 + Macro_Bio6 + Macro_Bio13 + Macro_Bio14 + landcover,
    data     = train_down,
    ntree    = 500,
    mtry     = 2,
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Lc_int_Macro <- predict(rf_Lc_int_Macro, test_data, type = "prob")[, "presence"]
  results$auc_rf_Lc_int_Macro[k] <- roc(test_data$pa, pred_rf_Lc_int_Macro, 
                                quiet = TRUE)$auc
  
  # ── 33. RF MacroCp * landcover ──────────────────────────────────────
  rf_Lc_int_MacroCp <- randomForest(
    pa_factor ~ Macro_CP_Tmin + Macro_CP_Tmax + Macro_CP_SMmax + Macro_CP_SMmin + landcover,
    data     = train_down,
    ntree    = 500,
    mtry     = 2,
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Lc_int_MacroCp <- predict(rf_Lc_int_MacroCp, test_data, type = "prob")[, "presence"]
  results$auc_rf_Lc_int_MacroCp[k] <- roc(test_data$pa, pred_rf_Lc_int_MacroCp, 
                                        quiet = TRUE)$auc
  
  # ── 34. RF MicroCp * landcover ──────────────────────────────────────
  rf_Lc_int_MicroCp <- randomForest(
    pa_factor ~ Micro_CP_Tmin + Micro_CP_Tmax + Micro_CP_SMmin + Micro_CP_SMmax + landcover,
    data     = train_down,
    ntree    = 500,
    mtry     = 2,
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Lc_int_MicroCp <- predict(rf_Lc_int_MicroCp, test_data, type = "prob")[, "presence"]
  results$auc_rf_Lc_int_MicroCp[k] <- roc(test_data$pa, pred_rf_Lc_int_MicroCp, 
                                          quiet = TRUE)$auc
  
  
  # ── 35. RF Microclimate * landcover ──────────────────────────────────────
  rf_Lc_int_Micro <- randomForest(
    pa_factor ~ Micro_Bio5 + Micro_Bio6 + Micro_Bio13 + Micro_Bio14 + landcover,
    data     = train_down,
    ntree    = 500,
    mtry     = 2,
    importance = TRUE,
    classwt    = c("background" = 1, "presence" = 1)
  )
  
  pred_rf_Lc_int_Micro <- predict(rf_Lc_int_Micro, test_data, type = "prob")[, "presence"]
  results$auc_rf_Lc_int_Micro[k] <- roc(test_data$pa, pred_rf_Lc_int_Micro, 
                                        quiet = TRUE)$auc
  
}

# Summarize the results, reporting the mean and sd
library(readxl)
Modeltable = read_xlsx("D:/Paper9/result/Models.xlsx")
summary_table <- data.frame(
  Model = Modeltable$`Model Name`,
  Mean_AUC = c(
    mean(results$auc_glm_Macro,       na.rm = TRUE),
    mean(results$auc_glm_MacroCp,       na.rm = TRUE),
    mean(results$auc_glm_MicroCp,   na.rm = TRUE),
    mean(results$auc_glm_Micro,   na.rm = TRUE),
    mean(results$auc_glm_Lc,   na.rm = TRUE),
    mean(results$auc_glm_Lc_int_Macro, na.rm = TRUE),
    mean(results$auc_glm_Lc_int_MacroCp, na.rm = TRUE),
    mean(results$auc_glm_Lc_int_MicroCp, na.rm = TRUE),
    mean(results$auc_glm_Lc_int_Micro, na.rm = TRUE),
    mean(results$auc_glm_Lc_add_Macro, na.rm = TRUE),
    mean(results$auc_glm_Lc_add_MacroCp, na.rm = TRUE),
    mean(results$auc_glm_Lc_add_MicroCp,   na.rm = TRUE),
    mean(results$auc_glm_Lc_add_Micro,   na.rm = TRUE),
    mean(results$auc_maxent_Macro, na.rm = TRUE),
    mean(results$auc_maxent_MacroCp, na.rm = TRUE),
    mean(results$auc_maxent_MicroCp, na.rm = TRUE),
    mean(results$auc_maxent_Micro, na.rm = TRUE),
    mean(results$auc_maxent_Lc, na.rm = TRUE),
    mean(results$auc_maxent_Lc_add_Macro, na.rm = TRUE),
    mean(results$auc_maxent_Lc_add_MacroCp, na.rm = TRUE),
    mean(results$auc_maxent_Lc_add_MicroCp, na.rm = TRUE),
    mean(results$auc_maxent_Lc_add_Micro, na.rm = TRUE),
    mean(results$auc_maxent_Lc_int_Macro, na.rm = TRUE),
    mean(results$auc_maxent_Lc_int_MacroCp, na.rm = TRUE),
    mean(results$auc_maxent_Lc_int_MicroCp, na.rm = TRUE),
    mean(results$auc_maxent_Lc_int_Micro,   na.rm = TRUE),
    mean(results$auc_rf_Macro,   na.rm = TRUE),
    mean(results$auc_rf_MacroCp, na.rm = TRUE),
    mean(results$auc_rf_MicroCp, na.rm = TRUE),
    mean(results$auc_rf_Micro, na.rm = TRUE),
    mean(results$auc_rf_Lc, na.rm = TRUE),
    mean(results$auc_rf_Lc_int_Macro, na.rm = TRUE),
    mean(results$auc_rf_Lc_int_MacroCp, na.rm = TRUE),
    mean(results$auc_rf_Lc_int_MicroCp, na.rm = TRUE),
    mean(results$auc_rf_Lc_int_Micro, na.rm = TRUE)
  ),
  SD_AUC = c(
    sd(results$auc_glm_Macro,       na.rm = TRUE),
    sd(results$auc_glm_MacroCp,       na.rm = TRUE),
    sd(results$auc_glm_MicroCp,   na.rm = TRUE),
    sd(results$auc_glm_Micro,   na.rm = TRUE),
    sd(results$auc_glm_Lc,   na.rm = TRUE),
    sd(results$auc_glm_Lc_int_Macro, na.rm = TRUE),
    sd(results$auc_glm_Lc_int_MacroCp, na.rm = TRUE),
    sd(results$auc_glm_Lc_int_MicroCp, na.rm = TRUE),
    sd(results$auc_glm_Lc_int_Micro, na.rm = TRUE),
    sd(results$auc_glm_Lc_add_Macro, na.rm = TRUE),
    sd(results$auc_glm_Lc_add_MacroCp, na.rm = TRUE),
    sd(results$auc_glm_Lc_add_MicroCp,   na.rm = TRUE),
    sd(results$auc_glm_Lc_add_Micro,   na.rm = TRUE),
    sd(results$auc_maxent_Macro, na.rm = TRUE),
    sd(results$auc_maxent_MacroCp, na.rm = TRUE),
    sd(results$auc_maxent_MicroCp, na.rm = TRUE),
    sd(results$auc_maxent_Micro, na.rm = TRUE),
    sd(results$auc_maxent_Lc, na.rm = TRUE),
    sd(results$auc_maxent_Lc_add_Macro, na.rm = TRUE),
    sd(results$auc_maxent_Lc_add_MacroCp, na.rm = TRUE),
    sd(results$auc_maxent_Lc_add_MicroCp, na.rm = TRUE),
    sd(results$auc_maxent_Lc_add_Micro, na.rm = TRUE),
    sd(results$auc_maxent_Lc_int_Macro, na.rm = TRUE),
    sd(results$auc_maxent_Lc_int_MacroCp, na.rm = TRUE),
    sd(results$auc_maxent_Lc_int_MicroCp, na.rm = TRUE),
    sd(results$auc_maxent_Lc_int_Micro,   na.rm = TRUE),
    sd(results$auc_rf_Macro,   na.rm = TRUE),
    sd(results$auc_rf_MacroCp, na.rm = TRUE),
    sd(results$auc_rf_MicroCp, na.rm = TRUE),
    sd(results$auc_rf_Micro, na.rm = TRUE),
    sd(results$auc_rf_Lc, na.rm = TRUE),
    sd(results$auc_rf_Lc_int_Macro, na.rm = TRUE),
    sd(results$auc_rf_Lc_int_MacroCp, na.rm = TRUE),
    sd(results$auc_rf_Lc_int_MicroCp, na.rm = TRUE),
    sd(results$auc_rf_Lc_int_Micro, na.rm = TRUE)
  )
)

summary_table$Mean_AUC <- round(summary_table$Mean_AUC, 4)
summary_table$SD_AUC   <- round(summary_table$SD_AUC,   4)
summary_table$Report   <- paste0(summary_table$Mean_AUC, 
                                 " ± ", summary_table$SD_AUC)

write.csv(summary_table, "../result/Model Output/Pitangus_sulphuratus_Auc.csv", row.names = FALSE)
write.csv(results, "../result/Model Output/Pitangus_sulphuratus_RawResults.csv")
