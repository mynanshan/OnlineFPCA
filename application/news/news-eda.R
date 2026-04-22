#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## ================================================================
## User settings
## ================================================================
project_root <- "."
data_dir <- file.path(project_root, "data", "news", "Data")
out_dir <- file.path(project_root, "data", "news", "processed")
out_file <- file.path(out_dir, "news_preprocessed.rds")
summary_file <- file.path(out_dir, "news_preprocess_summary.csv")
show_plots <- interactive()

## Public data provide three social platforms although Yang & Yao (2023)
## describe "four platforms" in the paper.
platforms <- c("Facebook", "GooglePlus", "LinkedIn")
threshold_total_clicks <- 6
morning_hour_lower <- 0
morning_hour_upper <- 6
clock_window <- c(6, 30)
minute_resolution <- 1L

## ================================================================
## Helpers
## ================================================================
find_metadata_file <- function(data_dir) {
  candidates <- c(
    file.path(data_dir, "News_Final.csv"),
    file.path(dirname(data_dir), "News_Final.csv"),
    file.path(data_dir, "news_final.csv"),
    file.path(dirname(data_dir), "news_final.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find News_Final.csv. Please check `data_dir`.")
  }
  hit[1]
}

find_feedback_files <- function(data_dir) {
  pats <- paste0("^(", paste(platforms, collapse = "|"), ")_.*\\.csv$")
  files <- list.files(data_dir, pattern = pats, full.names = TRUE)
  if (length(files) == 0) {
    stop("Cannot find any social feedback CSV files in `data_dir`.")
  }
  files
}

parse_publish_hour <- function(x, tz = "UTC") {
  tt <- as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = tz)
  as.numeric(format(tt, "%H")) +
    as.numeric(format(tt, "%M")) / 60 +
    as.numeric(format(tt, "%S")) / 3600
}

clock_label <- function(x) sprintf("%02d:00", as.integer(round(x)) %% 24)

standardize_cumulative <- function(ts_mat) {
  ts_mat <- as.matrix(ts_mat)
  storage.mode(ts_mat) <- "double"
  ts_mat[ts_mat < 0] <- 0
  t(apply(ts_mat, 1, cummax))
}

cumulative_to_hourly <- function(cum_144) {
  idx_hour <- seq(3, 144, by = 3)
  cum_48 <- cum_144[, idx_hour, drop = FALSE]
  hourly_48 <- cbind(cum_48[, 1, drop = FALSE], cum_48[, -1, drop = FALSE] - cum_48[, -ncol(cum_48), drop = FALSE])
  colnames(hourly_48) <- sprintf("H%02d", 1:48)
  hourly_48
}

read_feedback_file <- function(path) {
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  ts_cols <- sprintf("TS%d", 1:144)
  keep_cols <- c("IDLink", ts_cols)
  if (!all(keep_cols %in% names(dat))) {
    stop("Feedback file does not contain expected columns: ", basename(path))
  }
  dat <- dat[, keep_cols]
  dat$IDLink <- as.numeric(dat$IDLink)
  dat[, ts_cols] <- standardize_cumulative(dat[, ts_cols, drop = FALSE])
  dat
}

subset_dat <- function(dat, idx) {
  out <- lapply(dat, function(x) {
    if (is.list(x) && !is.data.frame(x)) x[idx] else x[idx, , drop = FALSE]
  })
  out
}

## ================================================================
## Read and aggregate raw data
## ================================================================
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

metadata_file <- find_metadata_file(data_dir)
feedback_files <- find_feedback_files(data_dir)

message("Reading metadata: ", metadata_file)
news_meta <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
news_meta$IDLink <- as.numeric(news_meta$IDLink)
news_meta$PublishDatePOSIX <- as.POSIXct(news_meta$PublishDate, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
news_meta$publish_hour <- parse_publish_hour(news_meta$PublishDate)

message("Reading and aggregating feedback trajectories from ", length(feedback_files), " files")
feedback_list <- lapply(feedback_files, read_feedback_file)
feedback_all <- bind_rows(feedback_list)
feedback_sum <- feedback_all %>%
  group_by(IDLink) %>%
  summarise(across(starts_with("TS"), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

hourly_48 <- cumulative_to_hourly(as.matrix(feedback_sum[, sprintf("TS%d", 1:144), drop = FALSE]))
feedback_hourly <- data.frame(IDLink = feedback_sum$IDLink, as.data.frame(hourly_48), check.names = FALSE)
feedback_sum$Total48 <- feedback_sum$TS144

news_joined <- news_meta %>%
  inner_join(feedback_sum %>% select(IDLink, Total48), by = "IDLink") %>%
  inner_join(feedback_hourly, by = "IDLink") %>%
  arrange(PublishDatePOSIX, IDLink)

news_joined <- news_joined %>%
  mutate(
    is_morning = publish_hour >= morning_hour_lower & publish_hour < morning_hour_upper,
    is_active = Total48 >= threshold_total_clicks
  )

morning_active <- news_joined %>%
  filter(is_morning, is_active) %>%
  arrange(PublishDatePOSIX, IDLink)

## ================================================================
## Build the online functional data objects
## ================================================================
minutes_total <- as.integer((clock_window[2] - clock_window[1]) * 60 / minute_resolution)
tgrid <- seq(0, 1, length.out = minutes_total + 1)

hour_cols <- sprintf("H%02d", 1:48)
hourly_mat <- as.matrix(morning_active[, hour_cols, drop = FALSE])

build_one_subject <- function(i) {
  pub_hour <- morning_active$publish_hour[i]
  obs_time_abs <- pub_hour + seq_len(48)
  keep <- obs_time_abs > clock_window[1] & obs_time_abs <= clock_window[2]

  y_hour <- hourly_mat[i, keep]
  t_abs <- obs_time_abs[keep]

  list(
    IDLink = morning_active$IDLink[i],
    PublishDate = morning_active$PublishDate[i],
    Topic = morning_active$Topic[i] %||% NA_character_,
    pub_hour = pub_hour,
    Total48 = morning_active$Total48[i],
    hourly24 = y_hour,
    time_abs = t_abs,
    Lt = (t_abs - clock_window[1]) / diff(clock_window),
    full_domain24 = length(y_hour) == 24
  )
}

subject_list <- lapply(seq_len(nrow(morning_active)), build_one_subject)
full_idx <- vapply(subject_list, `[[`, logical(1), "full_domain24")
subject_list <- subject_list[full_idx]

if (!all(full_idx)) {
  warning(sum(!full_idx), " morning-active rows did not yield 24 retained hourly measurements and were removed.")
}

hourly24_raw <- t(vapply(subject_list, `[[`, numeric(24), "hourly24"))
log24_raw <- log1p(hourly24_raw)
Lt_list <- lapply(subject_list, `[[`, "Lt")
time_abs_list <- lapply(subject_list, `[[`, "time_abs")

## Smooth mean for centering (raw log-click scale)
mean_fit <- smooth.spline(
  x = unlist(Lt_list),
  y = unlist(lapply(seq_len(nrow(log24_raw)), function(i) log24_raw[i, ])),
  nknots = min(100L, minutes_total + 1L)
)
mean_raw_eval <- predict(mean_fit, x = tgrid)$y

Ly_raw <- lapply(seq_len(nrow(log24_raw)), function(i) as.numeric(log24_raw[i, ]))
Ly_centered <- Map(function(y, tt) {
  y - predict(mean_fit, x = tt)$y
}, Ly_raw, Lt_list)

ysd <- sd(unlist(Ly_centered))
Ly_scaled <- lapply(Ly_centered, function(y) y / ysd)
Ltid_list <- lapply(Lt_list, function(tt) {
  pmax(1L, pmin(length(tgrid), round(tt * (length(tgrid) - 1)) + 1L))
})
Lmi_vec <- vapply(Ly_scaled, length, integer(1))

log24_centered <- t(vapply(Ly_centered, identity, numeric(24)))
log24_scaled <- t(vapply(Ly_scaled, identity, numeric(24)))

summary_tbl <- data.frame(
  step = c(
    "Rows in News_Final.csv",
    "Rows with aggregated feedback",
    "Morning news (0:00 <= publish < 6:00)",
    paste0("Morning + total clicks >= ", threshold_total_clicks),
    "Retained with 24 hourly observations on (6, 30]"
  ),
  count = c(
    nrow(news_meta),
    nrow(news_meta %>% inner_join(feedback_sum, by = "IDLink")),
    sum(news_joined$is_morning),
    nrow(morning_active),
    length(subject_list)
  )
)
write.csv(summary_tbl, summary_file, row.names = FALSE)

message("Preprocessing summary:")
print(summary_tbl)
message("Retained subjects: ", length(subject_list))
message("Paper target after filtering: 5689 morning-active news items.")

news_obj <- list(
  paper_reference = list(
    subset = "Morning news published between 0:00 and 6:00",
    clock_window = clock_window,
    dense_measurements = 24,
    threshold_total_clicks = threshold_total_clicks
  ),
  preprocessing = list(
    platforms = platforms,
    feedback_files = basename(feedback_files),
    metadata_file = basename(metadata_file),
    negative_value_handling = "Replace leading negative cumulative entries by 0, then enforce monotonicity with cummax.",
    hourly_aggregation = "Use TS3, TS6, ..., TS144 as cumulative hourly totals, then difference to obtain hourly clicks.",
    log_transform = "log1p(hourly_clicks)",
    mean_centering = "Global smoothing spline on pooled (time, log-click) pairs.",
    scaling = "Divide centered observations by pooled SD.",
    minute_resolution = minute_resolution
  ),
  summary = summary_tbl,
  t0 = 0,
  t1 = 1,
  tgrid = tgrid,
  clock_grid = clock_window[1] + diff(clock_window) * tgrid,
  mean_curve = data.frame(
    Lt = tgrid,
    clock = clock_window[1] + diff(clock_window) * tgrid,
    mean_raw = mean_raw_eval,
    mean_scaled = mean_raw_eval / ysd
  ),
  dat = list(
    IDLink = vapply(subject_list, `[[`, numeric(1), "IDLink"),
    PublishDate = vapply(subject_list, `[[`, character(1), "PublishDate"),
    Topic = vapply(subject_list, `[[`, character(1), "Topic"),
    pub_hour = vapply(subject_list, `[[`, numeric(1), "pub_hour"),
    Total48 = vapply(subject_list, `[[`, numeric(1), "Total48"),
    Ly = Ly_scaled,
    Ly_centered = Ly_centered,
    Ly_raw = Ly_raw,
    Ltid = Ltid_list,
    Lt = Lt_list,
    time_abs = time_abs_list,
    Lmi = as.integer(Lmi_vec)
  ),
  dense = list(
    hourly24_raw = hourly24_raw,
    log24_raw = log24_raw,
    log24_centered = log24_centered,
    log24_scaled = log24_scaled,
    clock_abs = do.call(rbind, time_abs_list),
    clock_rel = do.call(rbind, Lt_list)
  )
)

saveRDS(news_obj, out_file)
message("Saved preprocessed object to: ", out_file)

## ================================================================
## EDA for RStudio / interactive use
## ================================================================
plot_df_pub <- news_joined %>%
  mutate(status = case_when(
    is_morning & is_active ~ "retained before domain check",
    is_morning & !is_active ~ "morning but inactive",
    TRUE ~ "other publish times"
  ))

p_pubhour <- ggplot(plot_df_pub, aes(x = publish_hour, fill = status)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.7) +
  labs(
    x = "Publication hour",
    y = "Count",
    title = "Publication-time distribution of news items"
  ) +
  theme_minimal(base_size = 13)

p_totalclicks <- ggplot(morning_active, aes(x = log10(Total48 + 1))) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = log10(threshold_total_clicks + 1), linetype = 2, linewidth = 1) +
  labs(
    x = "log10(total clicks over 48h + 1)",
    y = "Count",
    title = "Activity filter for morning news"
  ) +
  theme_minimal(base_size = 13)

mean_plot_df <- news_obj$mean_curve
p_mean <- ggplot(mean_plot_df, aes(x = clock, y = mean_raw)) +
  geom_line(linewidth = 1.1) +
  scale_x_continuous(
    breaks = seq(clock_window[1], clock_window[2], by = 4),
    labels = clock_label
  ) +
  labs(
    x = "Clock time",
    y = "Estimated mean log-click count",
    title = "Mean curve on the Yang–Yao morning-news domain"
  ) +
  theme_minimal(base_size = 13)

set.seed(123)
nsamp <- min(120, nrow(log24_raw))
samp_idx <- sample(seq_len(nrow(log24_raw)), nsamp)
spaghetti_df <- data.frame(
  id = rep(seq_along(samp_idx), each = 24),
  clock = as.vector(do.call(rbind, time_abs_list[samp_idx])),
  log_click = as.vector(log24_raw[samp_idx, ])
)
p_spaghetti <- ggplot(spaghetti_df, aes(x = clock, y = log_click, group = id)) +
  geom_line(alpha = 0.18, linewidth = 0.4, color = "#1f77b4") +
  scale_x_continuous(
    breaks = seq(clock_window[1], clock_window[2], by = 4),
    labels = clock_label
  ) +
  labs(
    x = "Clock time",
    y = "log(hourly clicks + 1)",
    title = "Sample of retained morning-news trajectories"
  ) +
  theme_minimal(base_size = 13)

p_topics <- news_joined %>%
  count(Topic, sort = TRUE) %>%
  ggplot(aes(x = reorder(Topic, n), y = n)) +
  geom_col(fill = "grey65") +
  coord_flip() +
  labs(x = NULL, y = "Count", title = "News counts by topic") +
  theme_minimal(base_size = 13)

if (show_plots) {
  print(p_pubhour)
  print(p_totalclicks)
  print(p_mean)
  print(p_spaghetti)
  print(p_topics)
}

invisible(news_obj)
