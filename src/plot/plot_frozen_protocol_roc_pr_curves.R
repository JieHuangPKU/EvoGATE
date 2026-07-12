#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(cowplot)
  library(grid)
  library(gridExtra)
})

DEFAULT_BASE_DIR <- "/home/jiehuang/software/fungi/ProGATE_v2/outputs/frozen_protocol_benchmark_v2"
DEFAULT_OUTPUT_DIR <- "/home/jiehuang/software/fungi/ProGATE_v2/results/plots"
DEFAULT_FEATURE_SETTING <- "ORT_EXP_SUB"
DEFAULT_MODELS <- c("GraphSAGE", "N2V_MLP", "MLP", "NB", "SVM")
DEFAULT_SPECIES <- c("fgraminearum_newlabel", "scerevisiae", "celegans", "human", "dmelanogaster")
DEFAULT_LAYOUT_NROW <- 2L
DEFAULT_LAYOUT_NCOL <- 3L

MODEL_COLORS <- c(
  "GraphSAGE" = "#332288FF",
  "N2V_MLP" = "#117733FF",
  "MLP" = "#CC6677FF",
  "NB" = "#999933FF",
  "SVM" = "#88CCEEFF"
)

DISPLAY_NAME_MAP <- c(
  "fgraminearum_newlabel" = "F. graminearum",
  "fgraminearum_oldlabel" = "F. graminearum (old label)",
  "scerevisiae" = "S. cerevisiae",
  "celegans" = "C. elegans",
  "human" = "H. sapiens",
  "dmelanogaster" = "D. melanogaster"
)

REQUIRED_COLUMNS <- c(
  "canonical_gene_id", "split", "label", "is_labeled", "pred_score",
  "protocol", "species", "regime", "model", "feature_setting", "split_version"
)

parse_cli_args <- function() {
  option_list <- list(
    make_option("--base_dir", type = "character", default = DEFAULT_BASE_DIR),
    make_option("--feature_setting", type = "character", default = DEFAULT_FEATURE_SETTING),
    make_option("--models", type = "character", default = paste(DEFAULT_MODELS, collapse = ",")),
    make_option("--species", type = "character", default = paste(DEFAULT_SPECIES, collapse = ",")),
    make_option("--mode", type = "character", default = "pooled"),
    make_option("--split_version", type = "character", default = ""),
    make_option("--output_dir", type = "character", default = DEFAULT_OUTPUT_DIR),
    make_option("--layout_nrow", type = "integer", default = DEFAULT_LAYOUT_NROW),
    make_option("--layout_ncol", type = "integer", default = DEFAULT_LAYOUT_NCOL),
    make_option("--width", type = "double", default = 12),
    make_option("--height", type = "double", default = 7.2),
    make_option("--verbose", action = "store_true", default = FALSE)
  )
  parser <- OptionParser(option_list = option_list)
  args <- optparse::parse_args(parser)
  args$mode <- match.arg(args$mode, c("pooled", "mean_seed"))
  args$model_list <- trimws(strsplit(args$models, ",", fixed = TRUE)[[1]])
  args$species_list <- trimws(strsplit(args$species, ",", fixed = TRUE)[[1]])
  args$split_version <- trimws(args$split_version)
  args
}

log_message <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(..., "\n")
  }
}

parse_bool_string <- function(x) {
  normalized <- tolower(trimws(as.character(x)))
  normalized %in% c("true", "t", "1", "yes", "y")
}

derive_seed_from_path <- function(path) {
  matched <- regmatches(path, regexpr("run_[0-9]+", path))
  if (length(matched) == 0 || identical(matched, character(0)) || is.na(matched)) {
    return(NA_integer_)
  }
  as.integer(sub("run_", "", matched))
}

discover_prediction_files <- function(base_dir, verbose = FALSE) {
  if (!dir.exists(base_dir)) {
    stop(sprintf("Base directory does not exist: %s", base_dir), call. = FALSE)
  }
  files <- list.files(base_dir, pattern = "predictions.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("No predictions.tsv files found under %s", base_dir), call. = FALSE)
  }
  log_message(sprintf("[discover] found %d predictions.tsv files", length(files)), verbose = verbose)
  tibble::tibble(path = sort(files))
}

load_single_prediction_file <- function(path) {
  df <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE, col_types = readr::cols(.default = "c"))
  missing_columns <- setdiff(REQUIRED_COLUMNS, colnames(df))
  if (length(missing_columns) > 0) {
    stop(sprintf("Missing required columns in %s: %s", path, paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  df <- df %>%
    mutate(
      label = suppressWarnings(as.numeric(label)),
      pred_score = suppressWarnings(as.numeric(pred_score)),
      is_labeled_flag = parse_bool_string(is_labeled),
      seed = derive_seed_from_path(path),
      source_path = path,
      target_id = if_else(!is.na(protocol) & protocol != "", protocol, species)
    )
  if (any(is.na(df$pred_score))) {
    stop(sprintf("Non-numeric pred_score values found in %s", path), call. = FALSE)
  }
  df
}

load_predictions <- function(file_table, feature_setting, models, species, verbose = FALSE) {
  data_list <- lapply(file_table$path, load_single_prediction_file)
  predictions <- bind_rows(data_list) %>%
    filter(
      .data$split == "test",
      .data$is_labeled_flag,
      .data$feature_setting == .env$feature_setting,
      .data$model %in% .env$models,
      .data$target_id %in% .env$species
    ) %>%
    mutate(
      model = factor(.data$model, levels = .env$models),
      target_id = factor(.data$target_id, levels = .env$species)
    )
  if (nrow(predictions) == 0) {
    stop("No test + labeled prediction rows matched the requested feature_setting/models/species.", call. = FALSE)
  }
  if (any(is.na(predictions$label))) {
    bad_paths <- unique(predictions$source_path[is.na(predictions$label)])
    stop(sprintf("Non-numeric label values found after filtering in: %s", paste(bad_paths, collapse = ", ")), call. = FALSE)
  }
  if (any(is.na(predictions$pred_score))) {
    bad_paths <- unique(predictions$source_path[is.na(predictions$pred_score)])
    stop(sprintf("Non-numeric pred_score values found after filtering in: %s", paste(bad_paths, collapse = ", ")), call. = FALSE)
  }
  log_message(sprintf("[load] retained %d prediction rows", nrow(predictions)), verbose = verbose)
  predictions
}

validate_group_consistency <- function(predictions, split_version = "", verbose = FALSE) {
  split_table <- predictions %>%
    distinct(target_id, model, feature_setting, split_version)
  if (split_version != "") {
    filtered <- predictions %>% filter(.data$split_version == split_version)
    if (nrow(filtered) == 0) {
      stop(sprintf("No rows matched split_version=%s", split_version), call. = FALSE)
    }
    predictions <- filtered
  } else {
    bad_groups <- split_table %>%
      count(target_id, model, feature_setting, name = "n_split_versions") %>%
      filter(n_split_versions > 1)
    if (nrow(bad_groups) > 0) {
      details <- bad_groups %>%
        mutate(group = sprintf("%s / %s / %s", target_id, model, feature_setting)) %>%
        pull(group)
      stop(
        sprintf(
          "Found multiple split_version values within requested groups. Re-run with --split_version. Groups: %s",
          paste(details, collapse = "; ")
        ),
        call. = FALSE
      )
    }
  }
  chosen_split_versions <- predictions %>% distinct(target_id, model, feature_setting, split_version)
  log_message(sprintf("[validate] using split_version values: %s", paste(unique(chosen_split_versions$split_version), collapse = ", ")), verbose = verbose)
  predictions
}

roc_curve_points <- function(labels, scores) {
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  scores <- scores[ord]
  positives <- sum(labels == 1)
  negatives <- sum(labels == 0)
  if (positives == 0 || negatives == 0) {
    stop("ROC requires both positive and negative labels.", call. = FALSE)
  }
  tp <- cumsum(labels == 1)
  fp <- cumsum(labels == 0)
  tpr <- c(0, tp / positives, 1)
  fpr <- c(0, fp / negatives, 1)
  curve <- tibble::tibble(x = fpr, y = tpr) %>%
    group_by(x) %>%
    summarise(y = max(y), .groups = "drop") %>%
    arrange(x)
  auc <- sum(diff(curve$x) * (head(curve$y, -1) + tail(curve$y, -1)) / 2)
  list(curve = curve, metric = auc)
}

pr_curve_points <- function(labels, scores) {
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]
  positives <- sum(labels == 1)
  if (positives == 0) {
    stop("PR requires at least one positive label.", call. = FALSE)
  }
  tp <- cumsum(labels == 1)
  fp <- cumsum(labels == 0)
  recall <- tp / positives
  precision <- tp / pmax(tp + fp, 1)
  curve <- tibble::tibble(x = c(0, recall), y = c(1, precision)) %>%
    group_by(x) %>%
    summarise(y = max(y), .groups = "drop") %>%
    arrange(x)
  delta_recall <- diff(curve$x)
  auc <- sum(delta_recall * tail(curve$y, -1))
  list(curve = curve, metric = auc)
}

interpolate_curve <- function(curve_df, grid_values) {
  approx(
    x = curve_df$x,
    y = curve_df$y,
    xout = grid_values,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

compute_pooled_roc <- function(df) {
  roc_curve_points(df$label, df$pred_score)
}

compute_pooled_pr <- function(df) {
  pr_curve_points(df$label, df$pred_score)
}

compute_mean_seed_curve <- function(df, curve_type = c("roc", "pr"), grid_size = 201L) {
  curve_type <- match.arg(curve_type)
  seeds <- sort(unique(df$seed))
  if (length(seeds) == 0 || any(is.na(seeds))) {
    stop("mean_seed mode requires run directories with explicit numeric seeds.", call. = FALSE)
  }
  grid_values <- seq(0, 1, length.out = grid_size)
  curve_matrix <- lapply(seeds, function(seed_value) {
    seed_df <- df %>% filter(seed == seed_value)
    curve_obj <- if (curve_type == "roc") compute_pooled_roc(seed_df) else compute_pooled_pr(seed_df)
    tibble::tibble(
      seed = seed_value,
      x = grid_values,
      y = interpolate_curve(curve_obj$curve, grid_values),
      metric = curve_obj$metric
    )
  }) %>% bind_rows()
  summary_curve <- curve_matrix %>%
    group_by(x) %>%
    summarise(
      y = mean(y, na.rm = TRUE),
      y_sd = sd(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      y_sd = if_else(is.na(.data$y_sd), 0, .data$y_sd),
      y_lower = pmax(.data$y - .data$y_sd, 0),
      y_upper = pmin(.data$y + .data$y_sd, 1),
      y = pmin(pmax(.data$y, 0), 1)
    )
  list(
    curve = summary_curve,
    metric = mean(unique(curve_matrix$metric), na.rm = TRUE),
    metric_sd = sd(unique(curve_matrix$metric), na.rm = TRUE)
  )
}

prepare_plot_curve_data <- function(curves_df, species_id, curve_type, models) {
  raw_species_curves <- curves_df %>%
    filter(.data$target_id == species_id, .data$curve_type == curve_type) %>%
    mutate(model = factor(model, levels = models))

  duplicate_keys <- raw_species_curves %>%
    count(model, curve_type, x, name = "n_rows") %>%
    filter(n_rows > 1)
  if (nrow(duplicate_keys) > 0) {
    duplicate_summary <- duplicate_keys %>%
      mutate(key = sprintf("%s/%s/x=%.6f (n=%d)", model, curve_type, x, n_rows)) %>%
      pull(key)
    stop(
      sprintf(
        "Duplicate plotting x values detected for species=%s: %s",
        species_id,
        paste(duplicate_summary, collapse = "; ")
      ),
      call. = FALSE
    )
  }

  species_curves <- raw_species_curves %>%
    group_by(target_id, model, curve_type, x) %>%
    summarise(
      y = max(y, na.rm = TRUE),
      y_lower = if ("y_lower" %in% colnames(curves_df)) min(y_lower, na.rm = TRUE) else NA_real_,
      y_upper = if ("y_upper" %in% colnames(curves_df)) max(y_upper, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      y_lower = if_else(is.finite(.data$y_lower), .data$y_lower, NA_real_),
      y_upper = if_else(is.finite(.data$y_upper), .data$y_upper, NA_real_)
    ) %>%
    arrange(.data$model, .data$x)

  post_dedup_duplicates <- species_curves %>%
    count(model, curve_type, x, name = "n_rows") %>%
    filter(n_rows > 1)
  if (nrow(post_dedup_duplicates) > 0) {
    stop(
      sprintf("Plotting data still contains duplicate x values after deduplication for species=%s", species_id),
      call. = FALSE
    )
  }

  species_curves
}

compute_curves_for_species <- function(predictions, species_id, models, mode, verbose = FALSE) {
  species_df <- predictions %>% filter(.data$target_id == species_id)
  missing_models <- setdiff(models, unique(as.character(species_df$model)))
  if (length(missing_models) > 0) {
    log_message(sprintf("[species=%s] missing models for feature setting: %s", species_id, paste(missing_models, collapse = ", ")), verbose = verbose)
  }
  if (nrow(species_df) == 0) {
    return(list(curves = tibble::tibble(), stats = tibble::tibble(), missing_models = models, has_data = FALSE))
  }
  available_models <- intersect(models, unique(as.character(species_df$model)))
  curves <- list()
  stats <- list()
  for (model_name in available_models) {
    model_df <- species_df %>% filter(as.character(model) == model_name)
    if (mode == "pooled") {
      roc_obj <- compute_pooled_roc(model_df)
      pr_obj <- compute_pooled_pr(model_df)
      roc_curve_for_plot <- roc_obj$curve
      pr_curve_for_plot <- pr_obj$curve
    } else {
      roc_obj <- compute_mean_seed_curve(model_df, "roc")
      pr_obj <- compute_mean_seed_curve(model_df, "pr")
      roc_curve_for_plot <- roc_obj$curve
      pr_curve_for_plot <- pr_obj$curve
    }
    curves[[model_name]] <- list(
      roc = roc_curve_for_plot,
      pr = pr_curve_for_plot
    )
    stats[[model_name]] <- tibble::tibble(
      target_id = species_id,
      model = model_name,
      roc_metric = roc_obj$metric,
      pr_metric = pr_obj$metric
    )
  }
  curve_rows <- bind_rows(lapply(names(curves), function(model_name) {
    bind_rows(
      curves[[model_name]]$roc %>% mutate(curve_type = "roc", model = model_name),
      curves[[model_name]]$pr %>% mutate(curve_type = "pr", model = model_name)
    )
  }))
  list(curves = curve_rows, stats = bind_rows(stats), missing_models = missing_models, has_data = TRUE)
}

make_panel_title_label <- function(species_id) {
  display_name <- DISPLAY_NAME_MAP[[species_id]]
  if (is.null(display_name) || is.na(display_name)) {
    display_name <- species_id
  }
  display_name
}

make_metric_labels <- function(stats_df, models, metric_col, metric_prefix) {
  label_map <- setNames(models, models)
  for (model_name in intersect(models, stats_df$model)) {
      metric_value <- stats_df %>% filter(.data$model == model_name) %>% pull(.data[[metric_col]])
    if (length(metric_value) > 0) {
      label_map[[model_name]] <- sprintf("%s (%.3f)", model_name, metric_value[[1]])
    }
  }
  label_map
}

build_empty_panel <- function(species_id, x_label, y_label, x_limits, y_limits) {
  main_plot <- ggplot() +
    annotate("text", x = mean(x_limits), y = mean(y_limits), label = "No data", size = 4) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE) +
    labs(x = x_label, y = y_label, title = make_panel_title_label(species_id)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.title = element_text(face = "bold", size = 10),
      legend.position = "none"
    )
  main_plot
}

build_species_panel <- function(species_id, curves_df, stats_df, models, curve_type = c("roc", "pr"), mode) {
  curve_type <- match.arg(curve_type)
  species_curves <- prepare_plot_curve_data(curves_df, species_id, curve_type, models)
  species_stats <- stats_df %>% filter(.data$target_id == species_id)

  if (nrow(species_curves) == 0) {
    if (curve_type == "roc") {
      return(build_empty_panel(species_id, "False Positive Rate", "True Positive Rate", c(0, 1), c(0, 1)))
    }
    return(build_empty_panel(species_id, "Recall", "Precision", c(0, 1), c(0, 1)))
  }

  metric_col <- if (curve_type == "roc") "roc_metric" else "pr_metric"
  labels <- make_metric_labels(species_stats, models, metric_col, toupper(curve_type))
  available_models <- intersect(models, unique(as.character(species_curves$model)))
  color_values <- MODEL_COLORS[available_models]
  if (length(color_values) < length(available_models)) {
    fallback_colors <- scales::hue_pal()(length(available_models))
    names(fallback_colors) <- available_models
    color_values <- fallback_colors
  }

  x_label <- if (curve_type == "roc") "False Positive Rate" else "Recall"
  y_label <- if (curve_type == "roc") "True Positive Rate" else "Precision"
  p <- ggplot(species_curves, aes(x = x, y = y, color = model, group = model)) +
    geom_line(linewidth = 0.85, na.rm = TRUE) +
    scale_color_manual(values = color_values, breaks = available_models, labels = labels[available_models], drop = FALSE) +
    labs(x = x_label, y = y_label, color = NULL, title = make_panel_title_label(species_id)) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5, size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.title = element_text(face = "bold", size = 10, color = "black"),
      axis.text = element_text(size = 8, color = "black"),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.85), color = "grey70"),
      legend.key.size = unit(0.32, "cm"),
      legend.text = element_text(size = 7)
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)

  if (curve_type == "roc") {
    p <- p + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4)
  }

  p
}

adjust_layout <- function(species, layout_nrow, layout_ncol, verbose = FALSE) {
  capacity <- layout_nrow * layout_ncol
  if (length(species) <= capacity) {
    return(list(nrow = layout_nrow, ncol = layout_ncol))
  }
  new_nrow <- ceiling(length(species) / layout_ncol)
  log_message(
    sprintf("[layout] requested panels=%d exceed %dx%d; auto-adjusting nrow to %d", length(species), layout_nrow, layout_ncol, new_nrow),
    verbose = verbose
  )
  list(nrow = new_nrow, ncol = layout_ncol)
}

compose_multi_panel_plot <- function(panel_list, nrow, ncol) {
  total_slots <- nrow * ncol
  if (length(panel_list) < total_slots) {
    blank_panel <- ggplot() + theme_void()
    panel_list <- c(panel_list, replicate(total_slots - length(panel_list), blank_panel, simplify = FALSE))
  }
  if (length(panel_list) == 1) {
    return(panel_list[[1]])
  }
  gridExtra::arrangeGrob(grobs = panel_list, nrow = nrow, ncol = ncol)
}

build_output_stub <- function(args) {
  models_key <- if (identical(args$model_list, DEFAULT_MODELS)) "default_models" else paste(args$model_list, collapse = "-")
  species_key <- if (identical(args$species_list, DEFAULT_SPECIES)) "default_species" else paste(args$species_list, collapse = "-")
  sprintf("frozen_protocol_%s_%s_%s", args$feature_setting, args$mode, paste(c(models_key, species_key), collapse = "_"))
}

save_plot_bundle <- function(plot_obj, output_base, width, height, default_copy = NULL) {
  pdf_path <- sprintf("%s.pdf", output_base)
  png_path <- sprintf("%s.png", output_base)
  ggplot2::ggsave(pdf_path, plot = plot_obj, device = cairo_pdf, width = width, height = height, units = "in", dpi = 300)
  ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = 300)
  if (!is.null(default_copy)) {
    file.copy(pdf_path, sprintf("%s.pdf", default_copy), overwrite = TRUE)
    file.copy(png_path, sprintf("%s.png", default_copy), overwrite = TRUE)
  }
  list(pdf = normalizePath(pdf_path), png = normalizePath(png_path))
}

report_coverage <- function(predictions, requested_species, requested_models) {
  available <- predictions %>%
    distinct(target_id, model, split_version) %>%
    arrange(target_id, model)
  missing <- expand.grid(target_id = requested_species, model = requested_models, stringsAsFactors = FALSE) %>%
    anti_join(available, by = c("target_id", "model"))
  list(available = available, missing = missing)
}

main <- function() {
        args <- parse_cli_args()
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

  file_table <- discover_prediction_files(args$base_dir, verbose = args$verbose)
  predictions <- load_predictions(file_table, args$feature_setting, args$model_list, args$species_list, verbose = args$verbose)
  predictions <- validate_group_consistency(predictions, args$split_version, verbose = args$verbose)

  coverage <- report_coverage(predictions, args$species_list, args$model_list)
  if (nrow(coverage$missing) > 0) {
    warning(
      sprintf(
        "Missing species/model combinations for feature_setting=%s: %s",
        args$feature_setting,
        paste(sprintf("%s/%s", coverage$missing$target_id, coverage$missing$model), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  species_results <- lapply(args$species_list, function(species_id) {
    compute_curves_for_species(predictions, species_id, args$model_list, args$mode, verbose = args$verbose)
  })
  curves_df <- bind_rows(lapply(seq_along(species_results), function(idx) {
    species_results[[idx]]$curves %>% mutate(target_id = args$species_list[[idx]])
  }))
  stats_df <- bind_rows(lapply(species_results, `[[`, "stats"))

  layout <- adjust_layout(args$species_list, args$layout_nrow, args$layout_ncol, verbose = args$verbose)
  roc_panels <- lapply(args$species_list, function(species_id) {
    build_species_panel(species_id, curves_df, stats_df, args$model_list, curve_type = "roc", mode = args$mode)
  })
  pr_panels <- lapply(args$species_list, function(species_id) {
    build_species_panel(species_id, curves_df, stats_df, args$model_list, curve_type = "pr", mode = args$mode)
  })

  roc_plot <- compose_multi_panel_plot(roc_panels, layout$nrow, layout$ncol)
  pr_plot <- compose_multi_panel_plot(pr_panels, layout$nrow, layout$ncol)

  stub <- build_output_stub(args)
  is_default_run <- identical(args$feature_setting, DEFAULT_FEATURE_SETTING) &&
    identical(args$model_list, DEFAULT_MODELS) &&
    identical(args$species_list, DEFAULT_SPECIES) &&
    identical(args$mode, "pooled") &&
    identical(args$layout_nrow, DEFAULT_LAYOUT_NROW) &&
    identical(args$layout_ncol, DEFAULT_LAYOUT_NCOL)
  roc_paths <- save_plot_bundle(
    roc_plot,
    file.path(args$output_dir, sprintf("roc_curves_%s", stub)),
    width = args$width,
    height = args$height,
    default_copy = if (is_default_run) file.path(args$output_dir, "roc_curves_frozen_protocol_default") else NULL
  )
  pr_paths <- save_plot_bundle(
    pr_plot,
    file.path(args$output_dir, sprintf("pr_curves_%s", stub)),
    width = args$width,
    height = args$height,
    default_copy = if (is_default_run) file.path(args$output_dir, "pr_curves_frozen_protocol_default") else NULL
  )

  split_versions_used <- predictions %>%
    distinct(target_id, split_version) %>%
    arrange(target_id)
  cat(sprintf("ROC_PDF\t%s\n", roc_paths$pdf))
  cat(sprintf("ROC_PNG\t%s\n", roc_paths$png))
  cat(sprintf("PR_PDF\t%s\n", pr_paths$pdf))
  cat(sprintf("PR_PNG\t%s\n", pr_paths$png))
  cat("SPLIT_VERSIONS\n")
  print(split_versions_used)
  cat("AVAILABLE_COMBINATIONS\n")
  print(coverage$available)
  if (nrow(coverage$missing) > 0) {
    cat("MISSING_COMBINATIONS\n")
    print(coverage$missing)
  }
}

main()
