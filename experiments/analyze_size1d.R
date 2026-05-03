library(dplyr)
library(tidyr)
library(readr)

exprmt <- "size1d"
dirpath <- file.path("experiments", exprmt)


# RMSE table ----------------------------------------------------

res <- do.call(
  rbind,
  lapply(
    1:5, \(algo) {
      do.call(
        rbind,
        lapply(
          1:100, \(i) {
            seed <- i
            n_digit_seed <- ceiling(log10(seed + 1))
            seedtext <- paste0(paste(rep("0", 4 - n_digit_seed), collapse = ""), seed)
            filename <- paste0("simu_", exprmt, "_alg", algo, "_sd", seedtext, ".csv")
            tryCatch(
              {
                read_csv(file.path(dirpath, filename), show_col_types = FALSE)
              },
              error = function(e) NULL
            )
          }
        )
      )
    }
  )
)

res <- res |>
  mutate(
    seed = as.integer(seed),
    Ninit = as.integer(Ninit),
    npc = as.integer(npc),
    nBatch = as.integer(nBatch)
  ) |>
  mutate(
    RMSEphi1.avg = if_else(is.na(RMSEphi1.avg), RMSEphi1, RMSEphi1.avg),
    RMSEphi2.avg = if_else(is.na(RMSEphi2.avg), RMSEphi2, RMSEphi2.avg),
    RMSEphi3.avg = if_else(is.na(RMSEphi3.avg), RMSEphi3, RMSEphi3.avg)
  )

N <- 20000
Nsmalls <- c(200, 500, seq(1000, N, 1000))
meth_ord <- c(
  paste0("OnlineFPCA-", c("sgd", "adagrad")),
  "OnlineCov",
  "Batch-PACE", "Batch-FACE", "Batch-REML", "Batch-SOAP"
)

sumres <- res |>
  group_by(Method, noise, StepSize, Ninit, npc, initMethod, N) |>
  summarise_at(
    vars(Time, RMSEphi1.avg, RMSEphi2.avg, RMSEphi3.avg),
    list(mean = mean)
  ) |>
  ungroup() |>
  group_by(Method, noise, StepSize, Ninit, npc, initMethod) |>
  mutate(SumTime = cumsum(Time_mean)) |>
  ungroup() |>
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
      "OnlineCov" ~ "OnlineCov",
      "Batch-PACE" ~ "PACE",
      "Batch-FACE" ~ "FACE",
      "Batch-REML" ~ "REML",
      "Batch-SOAP" ~ "SOAP"
    ),
    SumTime = round(SumTime, 1)
  ) |>
  rename(
    RMSEphi1 = RMSEphi1.avg_mean,
    RMSEphi2 = RMSEphi2.avg_mean,
    RMSEphi3 = RMSEphi3.avg_mean
  )


if (interactive()) {
  View(
    res |> filter(N %in% Nsmalls) |>
      dplyr::select(
        noise, seed, Method, StepSize, N,
        tau, starts_with("RMSEphi")
      ),
    "Results"
  )
}


# plotting ----------------------------------------------------------------


library(ggplot2)
library(stringr)
library(scales)


plot_phi_rmse <- function(
    dat,
    noise_level = 0.5,
    x_var = c("N", "SumTime"),
    x_break = 2500,
    y_limits = NULL,
    show_initial_online = FALSE,
    save_path = NULL
) {

  x_var <- match.arg(x_var)

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
        labels = c("phi[1]", "phi[2]", "phi[3]")
      ),
      LegendLabel = as.character(Method)
    )

  if (!show_initial_online) {
    dat_long <- dat_long %>%
      filter(!(str_detect(LegendLabel, "^OnlineFPCA") & N == 0))
  }

  line_dat <- dat_long %>%
    filter(!is.na(RMSE)) %>%
    mutate(x_value = .data[[x_var]]) %>%
    arrange(component, LegendLabel, N)

  if (x_var == "SumTime") {
    line_dat <- line_dat %>%
      filter(x_value > 0)
  }

  x_end <- max(line_dat$x_value, na.rm = TRUE)
  if (is.null(y_limits)) {
    y_limits <- c(0, max(line_dat$RMSE, na.rm = TRUE) * 1.05)
  }

  if (x_var == "N") {
    x_limits <- c(0, x_end)
    x_breaks <- seq(0, x_end, by = x_break)
    x_labels <- label_number(scale = 1 / 1000)
    x_name <- "Number of processed subjects, N (x1000)"
    x_scale <- scale_x_continuous(
      name = x_name,
      limits = x_limits,
      breaks = x_breaks,
      labels = x_labels
    )
  } else {
    x_start <- min(line_dat$x_value, na.rm = TRUE)
    x_limits <- c(x_start, x_end)
    x_breaks <- scales::breaks_log(n = 6)(x_limits)
    x_breaks <- x_breaks[x_breaks >= x_start & x_breaks <= x_end]
    x_name <- "Cumulative computing time"
    x_scale <- scale_x_log10(
      name = x_name,
      limits = x_limits,
      breaks = x_breaks,
      labels = label_number()
    )
  }

  method_levels <- c(
    "OnlineFPCA-RSGD",
    "OnlineFPCA-RAdagrad",
    "OnlineCov",
    "PACE",
    "FACE",
    "REML",
    "SOAP"
  )

  color_values <- c(
    "OnlineFPCA-RSGD"     = "#0072B2",
    "OnlineFPCA-RAdagrad" = "#D55E00",
    "OnlineCov"           = "#009E73",
    "PACE"                = "#CC79A7",
    "FACE"                = "#E69F00",
    "REML"                = "#56B4E9",
    "SOAP"                = "black"
  )

  linetype_values <- c(
    "OnlineFPCA-RSGD"     = "solid",
    "OnlineFPCA-RAdagrad" = "longdash",
    "OnlineCov"           = "twodash",
    "PACE"                = "dashed",
    "FACE"                = "dotdash",
    "REML"                = "dotted",
    "SOAP"                = "solid"
  )

  active_method_levels <- method_levels[method_levels %in% unique(line_dat$LegendLabel)]

  p <- ggplot(line_dat) +
    geom_vline(
      xintercept = x_breaks[x_breaks > 0 & x_breaks <= x_end],
      linewidth = 0.25,
      color = "grey85"
    ) +

    geom_line(
      aes(
        x = x_value,
        y = RMSE,
        color = LegendLabel,
        linetype = LegendLabel
      ),
      linewidth = 0.85
    ) +

    scale_color_manual(
      values = color_values,
      breaks = active_method_levels,
      name = NULL
    ) +

    scale_linetype_manual(
      values = linetype_values,
      breaks = active_method_levels,
      name = NULL
    ) +
    x_scale +

    scale_y_continuous(
      name = expression("RMSE of estimated " * phi[k]),
      limits = y_limits,
      expand = expansion(mult = c(0.02, 0.08))
    ) +

    facet_wrap(
      ~ component,
      nrow = 1,
      labeller = label_parsed
    ) +

    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.title.position = "plot",

      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(size = 11),

      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),

      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.55, "cm"),

      axis.title.x.top = element_text(margin = margin(b = 6)),
      axis.title.x.bottom = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),

      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    ) +

    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linetype = unname(linetype_values[active_method_levels]),
          linewidth = rep(0.85, length(active_method_levels))
        )
      ),
      linetype = "none"
    )

  if (!is.null(save_path)) {
    ggsave(
      filename = save_path,
      plot = p,
      width = 8.4,
      height = 4.2,
      units = "in",
      device = cairo_pdf
    )
  }

  p
}

noise_sd <- 0.5
y_limits <- c(
  0,
  1.05 * max(unlist(sumres |> filter(noise == noise_sd) |> dplyr::select(starts_with("RMSEphi"))), na.rm = TRUE)
)
N_breaks <- 2500
p <- plot_phi_rmse(
  dat = sumres |>
  filter(N %% N_breaks == 0),
  noise_level = noise_sd,
  x_var = "N",
  x_break = N_breaks,
  y_limits = c(0, 0.5),
  save_path = NULL
)

p

ggsave(
  file.path(dirpath, "size1d_rmse.pdf"),
  p,
  width = 8.4,
  height = 4.2,
  dpi = 300,
  units = "in",
  device = cairo_pdf
)
