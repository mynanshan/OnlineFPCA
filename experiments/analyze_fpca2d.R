library(dplyr)
library(tidyr)
library(readr)

exprmt <- "fpca2d"
dirpath <- file.path("experiments", exprmt)


# RMSE table ----------------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:100, \(i) {
      seed <- i
      n_digit_seed <- ceiling(log10(seed + 1))
      seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
      filename <- paste0("simu_", exprmt, "_sd", seedtext, ".csv")
      tryCatch(
        read_csv(file.path(dirpath, filename), show_col_types = FALSE),
        error = function(e) NULL
      )
    }
  )
)

res <- res |>
  mutate(
    seed = as.integer(seed),
    nBatch = as.integer(nBatch)
  )

# Nsmall <- c(200, 300, 400, 500)
Nsmall <- c(200, 500)
meth_ord <- c("OnlineFPCA-sgd", "OnlineFPCA-adagrad", "Batch-SOAP", paste0("Batch-mOpCov-N", Nsmall))

sumres <- res |>
  group_by(Method, noise, StepSize, N) |>
  summarise_at(vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg), list(mean = mean)) |>
  ungroup() |>
  group_by(Method, noise, StepSize) |>
  mutate(SumTime = cumsum(Time_mean)) |> 
  filter(is.na(StepSize) | StepSize == 0.1) |>
  dplyr::select(
    Method, N, noise, StepSize, SumTime,
    RMSEphi1.avg_mean, RMSEphi2.avg_mean, RMSEphi3.avg_mean
  ) |>
  mutate(Method = factor(Method, levels = meth_ord)) |>
  relocate(noise, .before = Method) |>
  mutate(
    Method = recode_values(
      Method,
      "OnlineFPCA-sgd" ~ "OnlineFPCA-RSGD",
      "OnlineFPCA-adagrad" ~ "OnlineFPCA-RAdagrad",
      "Batch-SOAP" ~ "SOAP",
      "Batch-mOpCov-N200" ~ "mOpCov (N=200)",
      "Batch-mOpCov-N300" ~ "mOpCov (N=300)",
      "Batch-mOpCov-N400" ~ "mOpCov (N=400)",
      "Batch-mOpCov-N500" ~ "mOpCov (N=500)"
    ),
    SumTime = round(SumTime, 1)
  ) |>
  rename(
    RMSEphi1 = RMSEphi1.avg_mean,
    RMSEphi2 = RMSEphi2.avg_mean,
    RMSEphi3 = RMSEphi3.avg_mean
  )

q <- 5
N1 <- 5000
N0 <- N1 * (1:q)

View(
  res |> filter(N %in% N0) |>
    dplyr::select(noise, seed, Method, StepSize, N, starts_with("RMSEphi")),
  "Results"
)

tabres <- ungroup(sumres) |>
  filter(stringr::str_detect(Method, "Batch") | N %in% N0) |>
  arrange(noise, Method, N) |>
  mutate( epoch = as.integer(N / N1) ) |>
  relocate(epoch, .before = SumTime) |>
  dplyr::select(-N)

knitr::kable(
  tabres |>
    filter(is.na(StepSize) | StepSize == 0.1) |>
    filter(!stringr::str_detect(Method, "mOpCov")) |>
    dplyr::select(-StepSize),
  "latex",
  digits = 3
)



# plotting ----------------------------------------------------------------


library(tidyverse)
library(scales)


plot_phi_rmse <- function(
    dat,
    noise_level = 0.1,
    n_subjects = 5000,
    show_initial_online = FALSE,
    save_path = NULL
) {
  
  dat_long <- dat %>%
    filter(near(noise, noise_level)) %>%
    pivot_longer(
      cols = starts_with("RMSEphi"),
      names_to = "component",
      values_to = "RMSE"
    ) %>%
    mutate(
      component_id = as.integer(str_extract(component, "\\d+")),
      component = factor(
        component_id,
        levels = 1:3,
        labels = c("\u03D5\u2081", "\u03D5\u2082", "\u03D5\u2083")
      ),
      
      method_type = case_when(
        str_detect(Method, "OnlineFPCA") ~ "online",
        TRUE ~ "batch"
      ),
      
      MethodLabel = case_when(
        Method == "SOAP" ~ paste0("SOAP-2D (N=", N, ")"),
        str_detect(Method, "^mOpCov") ~ paste0("mOpCov (N=", N, ")"),
        TRUE ~ Method
      ),
      
      LegendLabel = case_when(
        Method == "SOAP" ~ "SOAP-2D",
        str_detect(Method, "^mOpCov") ~ "mOpCov",
        TRUE ~ Method
      )
    )
  
  if (!show_initial_online) {
    dat_long <- dat_long %>%
      filter(!(method_type == "online" & N == 0))
  }
  
  online_dat <- dat_long %>%
    filter(method_type == "online", !is.na(RMSE)) %>%
    arrange(component, MethodLabel, N)
  
  batch_dat <- dat_long %>%
    filter(method_type == "batch", !is.na(RMSE)) %>%
    mutate(
      N_label = str_extract(MethodLabel, "N=\\d+"),
      point_shape = if_else(str_detect(MethodLabel, "SOAP"), 15L, 21L)
    )
  
  x_end <- max(dat_long$N, na.rm = TRUE)
  
  batch_dat <- batch_dat %>%
    mutate(xend = x_end)
  
  method_levels <- c(
    "OnlineFPCA-RSGD",
    "OnlineFPCA-RAdagrad",
    "SOAP-2D",
    "mOpCov"
  )
  
  color_values <- c(
    "OnlineFPCA-RSGD"     = "#0072B2",
    "OnlineFPCA-RAdagrad" = "#D55E00",
    "SOAP-2D"             = "black",
    "mOpCov"              = "grey35"
  )
  
  linetype_values <- c(
    "OnlineFPCA-RSGD"     = "solid",
    "OnlineFPCA-RAdagrad" = "longdash",
    "SOAP-2D"             = "dashed",
    "mOpCov"              = "dotted"
  )
  
  shape_values <- c(
    "OnlineFPCA-RSGD"     = NA,
    "OnlineFPCA-RAdagrad" = NA,
    "SOAP-2D"             = 15,
    "mOpCov"              = 21
  )
  
  p <- ggplot() +
    geom_vline(
      xintercept = seq(n_subjects, x_end, by = n_subjects),
      linewidth = 0.25,
      color = "grey85"
    ) +
    
    geom_segment(
      data = batch_dat,
      aes(
        x = N, xend = xend,
        y = RMSE, yend = RMSE,
        color = LegendLabel,
        linetype = LegendLabel
      ),
      linewidth = 0.45,
      alpha = 0.85
    ) +
    
    geom_line(
      data = online_dat,
      aes(
        x = N,
        y = RMSE,
        color = LegendLabel,
        linetype = LegendLabel
      ),
      linewidth = 0.85
    ) +
    
    geom_point(
      data = batch_dat,
      aes(
        x = N,
        y = RMSE,
        color = LegendLabel,
        # shape = point_shape
        shape = LegendLabel
      ),
      size = 2.3,
      stroke = 0.8,
      fill = "white",
      show.legend = TRUE
    ) +
    
    # scale_shape_identity() +
    scale_shape_manual(
      values = shape_values,
      breaks = method_levels,
      name = NULL,
      drop = FALSE
    ) +
    
    scale_color_manual(
      values = color_values,
      breaks = method_levels,
      name = NULL
    ) +
    
    scale_linetype_manual(
      values = linetype_values,
      breaks = method_levels,
      name = NULL
    ) +
    
    scale_x_continuous(
      name = "Number of processed subjects, N",
      limits = c(0, x_end),
      breaks = seq(0, x_end, by = n_subjects),
      labels = label_comma(),
      sec.axis = sec_axis(
        trans = ~ . / n_subjects,
        breaks = seq(0, x_end / n_subjects, by = 1),
        name = "Epoch"
      )
    ) +
    
    scale_y_continuous(
      name = "RMSE of estimated \u03D5\u2096",
      limits = c(0, 0.25),
      expand = expansion(mult = c(0.02, 0.08))
    ) +
    
    facet_wrap(
      ~ component,
      ncol = 1
    ) +
    
    labs(
      title = paste0("Noise level: ", noise_level)
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.title.position = "plot",
      
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(size = 11),
      
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.55, "cm"),
      
      axis.title.x.top = element_text(margin = margin(b = 6)),
      axis.title.x.bottom = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    
    # guides(
    #   color = guide_legend(
    #     ncol = 1,
    #     byrow = TRUE,
    #     override.aes = list(
    #       linetype = c("solid", "longdash", "dashed", "dotted"),
    #       shape    = c(NA, NA, 15, 21),
    #       fill     = c(NA, NA, "black", "white"),
    #       linewidth = c(0.85, 0.85, 0.45, 0.45),
    #       size = c(NA, NA, 2.3, 2.3)
    #     )
    #   ),
    #   linetype = "none"
    # )
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linetype  = unname(linetype_values[method_levels]),
          shape     = unname(shape_values[method_levels]),
          linewidth = c(0.85, 0.85, 0.45, 0.45),
          size      = c(NA, NA, 2.3, 2.3),
          fill      = c(NA, NA, "black", "white"),
          alpha     = 1
        )
      ),
      linetype = "none",
      shape = "none"
    )
  
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = batch_dat,
        aes(
          x = N,
          y = RMSE,
          label = N_label,
          color = LegendLabel
        ),
        size = 2.6,
        min.segment.length = 0,
        box.padding = 0.15,
        point.padding = 0.15,
        nudge_x = 0.015 * x_end,
        direction = "y",
        show.legend = FALSE,
        max.overlaps = Inf
      )
  } else {
    p <- p +
      geom_text(
        data = batch_dat,
        aes(
          x = N,
          y = RMSE,
          label = N_label,
          color = LegendLabel
        ),
        size = 2.6,
        hjust = -0.1,
        vjust = -0.4,
        show.legend = FALSE
      )
  }
  
  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = p,
      width = 6.8,
      height = 8.4,
      units = "in",
      device = cairo_pdf
    )
  }
  
  return(p)
}

p1 <- plot_phi_rmse(dat = sumres, noise_level = 0.1, save_path = NULL)
p2 <- plot_phi_rmse(dat = sumres, noise_level = 0.5, save_path = NULL)
p3 <- plot_phi_rmse(dat = sumres, noise_level = 1.0, save_path = NULL)

library(patchwork)

p <- (p1 + p2 + p3) +
  plot_layout(
    ncol = 3,
    guides = "collect",
    axes = "collect_y",
    axis_titles = "collect_y",
    widths = c(1, 1, 1)
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.justification = "center",
    legend.key.width = unit(1.4, "cm")
  )

p

ggsave(
  file.path(dirpath, "fpca2d_rmse.pdf"),
  p,
  width = 8.4,
  height = 9,
  dpi = 300,
  units = "in",
  device = cairo_pdf
)
