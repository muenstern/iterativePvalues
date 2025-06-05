#' Prepare data (against zero)

#' Prepares dataframes for testing against zero
#' @param df_input input dataframe
#' @param var dependent variable
#' @param id unique identifier
#' @return Prepared dataframe
#' @export

prepare_data_zero <- function(df_input, var, id) {
  df_prepared <- data.frame()
  df_input$dummy <- var
  df_input$zero <- 0
  df_input <- df_input %>% rename("dv" = var)

  df_input$dummy_zero <- paste0(var, "_zero")

  df1 <- df_input[, c("VP", "dummy", "dv")]
  colnames(df1) <- c("VP", "dummy", "dv")

  df2 <- df_input[, c("VP", "dummy_zero", "zero")]
  colnames(df2) <- c("VP", "dummy", "dv")

  df_prepared <<- rbind(df1, df2)
  df_prepared <- as.data.frame(df_prepared)
  list2env(df_prepared, envir = .GlobalEnv)
  df_prepared <- as.data.frame(df_prepared)
  return(df_prepared)
}

#' Prepare data (against each other)

#' #' Prepares dataframes for testing against each other
#' @param df_input input dataframe
#' @param factor two-staged factor for testing against each other
#' @param var dependent variable
#' @param id unique identifier
#' @return prepared dataframe
#' @export

prepare_data_against <- function(df_input, factor, var, id) {
  df_prepared <- data.frame()
  df_input <- df_input %>% ungroup %>% select(id, factor, var)
  df_input <- df_input %>% rename("dummy" = factor)
  df_input <- df_input %>% rename("dv" = var)
  df_prepared <<- df_input
  df_prepared <- as.data.frame(df_prepared)
  list2env(df_prepared, envir = .GlobalEnv)
  df_prepared <- as.data.frame(df_prepared)
  return(df_prepared)
}


#' Iterative p-values

#' Computes iterative p-values
#' @param df input dataframe from prepare data function
#' @param n_iter number of iterations (default value = 100)
#' @param fct_levels default value = 1
#' @return dataframe "all_iterations"
#' @export

compute_p_curves <- function(n_iter = 100, df = df_prepared, fct_levels = 1) {
  n_iterations <<- n_iter
  all_iterations <- data.frame()
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  for (iteration in 1:n_iter) {
    setTxtProgressBar(pb, iteration)

    unique_VPs <- sample(unique(df$VP))
    get_appended <- data.frame()
    appended_pvalues <- data.frame()
    n <- 0

    for (i in unique_VPs) {
      if (i %in% get_appended$VP) next

      df_single_VP <- df %>% filter(VP == i)
      get_appended <- bind_rows(get_appended, df_single_VP)

      valid_grouping <- get_appended %>%
        group_by(dummy) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n > 1)

      if (nrow(valid_grouping) < 2 || length(unique(get_appended$VP)) < fct_levels) {
        p <- NA
      } else {
        t_test_result <- tryCatch(
          t.test(dv ~ dummy, data = get_appended, paired = TRUE),
          error = function(e) NULL
        )
        p <- ifelse(is.null(t_test_result), NA, t_test_result$p.value)
      }

      n <- n + 1

      appended_pvalues <- bind_rows(appended_pvalues, data.frame(
        N = n,
        VP = i,
        p_value = p,
        iteration = iteration
      ))
    }

    all_iterations <- bind_rows(all_iterations, appended_pvalues)
  }

  all_iterations <<- all_iterations

  close(pb)

  list2env(all_iterations, envir = .GlobalEnv)
  return(all_iterations)
}

#' Summary

#' Summarizes information about iterative p-value-function
#' @param df_iterations input dataframe from compute_p_curves function
#' @return dataframe "summary_stats"
#' @export

summarize_p_values <- function(df_iterations = all_iterations) {
  summary_stats <- data.frame()
  summary_stats <<- all_iterations %>%
    group_by(N) %>%
    summarise(
      mean_p = mean(p_value, na.rm = TRUE),
      sd_p = sd(p_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      min_p_value = mean_p - sd_p,
      max_p_value = mean_p + sd_p,
      abs_diff_sd_p = abs(c(NA, diff(sd_p)))
    ) %>%
    drop_na()

  return(summary_stats)
  list2env(summary_stats, envir = .GlobalEnv)
}


#' Parameters and Plot

#' Gives iterative-p-value plot and metrics
#' @param summarized_data data from summarize_p_values-function
#' @param n_iterations number of iterations
#' @return parameters and plot
#' @export

p_plot <- function(summarized_data = summary_stats){

  n_iter <- n_iterations

  # Regressionsanalysen (ungeglättete Werte) #
  lm_model <- lm(mean_p ~ N, data = summarized_data)
  summary_lm <- summary(lm_model)

  intercept <- coef(lm_model)[1]
  slope <<- coef(lm_model)[2]
  r_squared <<- summary_lm$r.squared
  p_val_slope <<- summary_lm$coefficients["N", "Pr(>|t|)"]
  beta_std <<- coef(lm(scale(mean_p) ~ scale(N), data = summarized_data))[2]

  # Metriken berechnen #
  auc_value <<- pracma::trapz(summarized_data$N, summarized_data$mean_p)
  tts_value <<- summarized_data$N[which(summarized_data$mean_p < 0.05)[1]]
  tts_max_value <<- summarized_data$N[which(summarized_data$max_p_value < 0.05)[1]]

  n_Total <<- length(unique(df_prepared$VP))
  sd_sig_level <- 1 - tts_max_value/n_Total
  sd_sig_level <<- round(sd_sig_level, 3) * 100

  auc_var <<- pracma::trapz(summarized_data$N, summarized_data$sd_p)
  auc_grad_var <<- pracma::trapz(summarized_data$N, summarized_data$abs_diff_sd_p)

  max_N <- max(summarized_data$N)

  # Normierte AUCs
  auc_norm <<- auc_value / max_N
  auc_var_norm <<- auc_var / max_N
  auc_grad_var_norm <<- auc_grad_var / max_N

  # Robuste Normierung über 95%-Quantile
  max_mean_p <- quantile(summarized_data$mean_p, 0.95, na.rm = TRUE)
  max_sd_p <- quantile(summarized_data$sd_p, 0.95, na.rm = TRUE)
  max_abs_diff_sd_p <- quantile(summarized_data$abs_diff_sd_p, 0.95, na.rm = TRUE)

  auc_double_norm <- pracma::trapz(summarized_data$N, summarized_data$mean_p / max_mean_p) / max_N
  auc_var_double_norm <- pracma::trapz(summarized_data$N, summarized_data$sd_p / max_sd_p) / max_N
  auc_grad_var_double_norm <- pracma::trapz(summarized_data$N, summarized_data$abs_diff_sd_p / max_abs_diff_sd_p) / max_N


  # Teilnehmerzahl extrahieren für Plot-Untertitel
  n_Total <- length(unique(df_prepared$VP))

  metrics <<- data.frame(TTS_M = numeric(), TTS_SD = numeric(), SD_05 = numeric(), AUC_p = numeric(), AUC_var = numeric(), AUC_gradvar = numeric())
  metrics <<- data.frame(TTS_M = tts_value, TTS_SD = tts_max_value, SD_05 = sd_sig_level, AUC_p = auc_value, AUC_var = auc_var, AUC_gradvar = auc_grad_var)

  p_plot <- ggplot(summarized_data, aes(x = N)) +
    theme_minimal(base_family = "Aptos") +

    geom_line(aes(y = mean_p, color = "Mittlerer p-Wert"), size = 1.2) +
    geom_ribbon(aes(ymin = min_p_value, ymax = max_p_value), fill = "#0072B2", alpha = 0.2) +
    geom_abline(intercept = intercept, slope = slope, linetype = "dotted", color = "black") +

    geom_hline(yintercept = 0.05, linetype = "solid", color = "red") +
    geom_vline(xintercept = tts_value, linetype = "dashed", color = "darkgreen") +
    geom_vline(xintercept = tts_max_value, linetype = "dashed", color = "darkgreen") +

    geom_line(aes(y = sd_p * 3, color = "SD der p-Werte"), linetype = "longdash", size = 0.8) +

    # TTS-Annotation
    annotate("text", x = tts_value, y = max(summarized_data$max_p_value), hjust = 1, vjust = 1, size = 5, family = "Cambria",
             label = paste0("TTS[M] = ", tts_value)) +
    annotate("text", x = tts_max_value, y = max(summarized_data$max_p_value)*0.8, hjust = 1, vjust = 1, size = 5, family = "Cambria",
             label = paste0("TTS[SD] = ", tts_max_value)) +

    scale_y_continuous(
      name = "p-values (M) ± SD",
      limits = c(-0.1, 1),
      breaks = seq(0, 1, 0.05),
      sec.axis = sec_axis(~./3, name = "p-values (SD)")
    ) +

    scale_color_manual(name = NULL,
                       values = c("p-values (M)" = "#0072B2", "p-values (SD)" = "darkorange")) +

    labs(
      title = paste0("Mean p-values over ", n_iter, " Iterationen"),
      subtitle = paste0(n_Total, " steps"),
      x = "Steps (n)"
    ) +

    theme(
      text = element_text(size = 20),
      legend.position = "bottom",
      legend.text = element_text(size = 14)
    )


  print(p_plot)

  cat("\n--- p-Value Plot Parameters ---\n")
  cat(sprintf("Iterationen = %s\n", n_iter))
  cat(sprintf("VP = %s\n", n_Total))

  cat("\n--- Regression Metrics ---\n")
  cat(sprintf("p = %.3f %s %.4f * N\n",
              intercept,
              ifelse(slope < 0, "-", "+"),
              abs(slope)))
  cat(sprintf("R² = %.3f, β = %.3f, p[slope] = %s\n\n",
              r_squared,
              beta_std,
              ifelse(p_val_slope < 0.001, "< .001", round(p_val_slope, 3))))

  cat("--- TTS ---\n")
  cat(sprintf("TTS[M] = %s\n", tts_value))
  cat(sprintf("TTS[SD] = %s\n", tts_max_value))
  cat(sprintf("SD < .05: %s / %s (%.2f%%)\n\n",
              tts_max_value,
              n_Total,
              sd_sig_level * 100))

  cat("--- AUC (absolute) ---\n")
  cat(sprintf("AUC[p]       = %.2f\n", auc_value))
  cat(sprintf("AUC[Var]     = %.4f\n", auc_var))
  cat(sprintf("AUC[GradVar] = %.4f\n\n", auc_grad_var))

  cat("--- AUC (standardized (N) ---\n")
  cat(sprintf("AUC[p]       = %.3f\n", auc_norm))
  cat(sprintf("AUC[Var]     = %.3f\n", auc_var_norm))
  cat(sprintf("AUC[GradVar] = %.3f\n", auc_grad_var_norm))

}

