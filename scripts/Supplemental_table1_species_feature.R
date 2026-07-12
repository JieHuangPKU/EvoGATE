#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(readr)
  library(yaml)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_path <- file.path(project_root, "configs", "frozen_protocol.yaml")
output_dir <- file.path(project_root, "results", "Supplemental")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- yaml::read_yaml(config_path)

object_order <- c(
  "fgraminearum_newlabel",
  "fgraminearum_oldlabel",
  "scerevisiae",
  "human",
  "celegans",
  "dmelanogaster"
)
object_labels <- c(
  fgraminearum_newlabel = "F. graminearum (new label)",
  fgraminearum_oldlabel = "F. graminearum (old label)",
  scerevisiae = "S. cerevisiae",
  human = "H. sapiens",
  celegans = "C. elegans",
  dmelanogaster = "D. melanogaster"
)
feature_rows <- c(
  orthologs = "Ortholog",
  expression = "Expression",
  sublocalization = "Localization"
)

format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

normalize_rel <- function(path) {
  sub(
    paste0("^", normalizePath(project_root, winslash = "/", mustWork = TRUE), "/?"),
    "",
    normalizePath(path, winslash = "/", mustWork = TRUE)
  )
}

read_header_width <- function(path) {
  header <- names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
  if (length(header) < 2) {
    stop(sprintf("Expected at least 2 columns in %s", path))
  }
  length(header) - 1L
}

read_gene_column <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  unique(as.character(df[[1]]))
}

read_ppi_stats <- function(path) {
  ppi <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (!all(c("A", "B") %in% names(ppi))) {
    stop(sprintf("Missing A/B columns in %s", path))
  }
  list(
    edge_count = nrow(ppi),
    ppi_node_ids = unique(c(as.character(ppi$A), as.character(ppi$B)))
  )
}

read_label_split_stats <- function(label_path, split_path) {
  labels <- readr::read_tsv(label_path, show_col_types = FALSE, progress = FALSE)
  splits <- readr::read_tsv(split_path, show_col_types = FALSE, progress = FALSE)
  list(
    labeled_count = nrow(labels),
    essential_count = sum(suppressWarnings(as.integer(labels$label)) == 1L, na.rm = TRUE),
    train_count = sum(splits$split == "train", na.rm = TRUE),
    val_count = sum(splits$split == "val", na.rm = TRUE),
    test_count = sum(splits$split == "test", na.rm = TRUE)
  )
}

read_object_info <- function(object_id) {
  protocol_cfg <- config$protocols[[object_id]]
  if (is.null(protocol_cfg)) {
    stop(sprintf("Protocol %s not found in %s", object_id, config_path))
  }

  data_key <- protocol_cfg$data_key
  ortholog_path <- file.path(project_root, config$feature_roots$orthologs, data_key, "orthologs.csv")
  expression_path <- file.path(project_root, config$feature_roots$expression, data_key, "profile.csv")
  subloc_path <- file.path(project_root, config$feature_roots$sublocalization, data_key, "subloc.csv")
  ppi_path <- file.path(project_root, config$feature_roots$ppi, data_key, "string.csv")
  label_path <- file.path(project_root, config$paths$labels_dir, protocol_cfg$label_output)
  split_path <- file.path(project_root, config$paths$splits_dir, protocol_cfg$split_output)

  required_paths <- c(ortholog_path, expression_path, subloc_path, ppi_path, label_path, split_path)
  missing <- required_paths[!file.exists(required_paths)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required files for %s:\n%s", object_id, paste(missing, collapse = "\n")))
  }

  feature_dimensions <- c(
    orthologs = read_header_width(ortholog_path),
    expression = read_header_width(expression_path),
    sublocalization = read_header_width(subloc_path)
  )
  ortholog_genes <- read_gene_column(ortholog_path)
  expression_genes <- read_gene_column(expression_path)
  subloc_genes <- read_gene_column(subloc_path)
  ppi_stats <- read_ppi_stats(ppi_path)
  label_stats <- read_label_split_stats(label_path, split_path)

  list(
    object_id = object_id,
    display_name = unname(object_labels[[object_id]]),
    regime = protocol_cfg$regime,
    species = protocol_cfg$species,
    data_key = data_key,
    feature_dimensions = feature_dimensions,
    node_count = length(unique(c(ortholog_genes, expression_genes, subloc_genes, ppi_stats$ppi_node_ids))),
    edge_count = ppi_stats$edge_count,
    labeled_count = label_stats$labeled_count,
    essential_count = label_stats$essential_count,
    train_count = label_stats$train_count,
    val_count = label_stats$val_count,
    test_count = label_stats$test_count,
    files = list(
      orthologs = normalize_rel(ortholog_path),
      expression = normalize_rel(expression_path),
      sublocalization = normalize_rel(subloc_path),
      ppi = normalize_rel(ppi_path),
      labels = normalize_rel(label_path),
      splits = normalize_rel(split_path)
    )
  )
}

dataset_info <- lapply(object_order, read_object_info)
names(dataset_info) <- object_order

feature_table <- data.frame("Feature Type" = unname(feature_rows), check.names = FALSE)
for (object_id in object_order) {
  feature_table[[object_labels[[object_id]]]] <- unname(dataset_info[[object_id]]$feature_dimensions[names(feature_rows)])
}

network_table <- data.frame(
  Statistic = c("N.nodes", "N.edges", "N.labeled genes", "N.essential genes", "N.train genes", "N.val genes", "N.test genes"),
  check.names = FALSE
)
for (object_id in object_order) {
  info <- dataset_info[[object_id]]
  network_table[[object_labels[[object_id]]]] <- c(
    info$node_count,
    info$edge_count,
    info$labeled_count,
    info$essential_count,
    info$train_count,
    info$val_count,
    info$test_count
  )
}

write_excel <- function(feature_df, network_df, output_path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Feature Dimensions", gridLines = FALSE)
  openxlsx::addWorksheet(wb, "Network Statistics", gridLines = FALSE)

  title_style_blue <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold", halign = "left")
  title_style_red <- openxlsx::createStyle(fontSize = 14, textDecoration = "bold", halign = "left")
  header_style_blue <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#4E79A7", textDecoration = "bold",
    halign = "center", border = "Bottom"
  )
  header_style_red <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#E15759", textDecoration = "bold",
    halign = "center", border = "Bottom"
  )
  first_col_style <- openxlsx::createStyle(textDecoration = "bold")
  striped_fill <- c("#FFFFFF", "#F5F7FA")

  write_sheet <- function(sheet_name, title, df, header_style, title_style) {
    openxlsx::writeData(wb, sheet_name, title, startRow = 1, startCol = 1, colNames = FALSE)
    openxlsx::addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE)
    openxlsx::writeData(wb, sheet_name, df, startRow = 3, startCol = 1, headerStyle = header_style)
    openxlsx::addStyle(wb, sheet_name, first_col_style, rows = 4:(nrow(df) + 3), cols = 1, gridExpand = TRUE, stack = TRUE)
    for (i in seq_len(nrow(df))) {
      row_index <- i + 3
      fill_style <- openxlsx::createStyle(fgFill = striped_fill[((i - 1) %% 2) + 1])
      openxlsx::addStyle(wb, sheet_name, fill_style, rows = row_index, cols = seq_len(ncol(df)), gridExpand = TRUE, stack = TRUE)
    }
    openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    openxlsx::freezePane(wb, sheet_name, firstActiveRow = 4)
  }

  write_sheet(
    "Feature Dimensions",
    "Supplemental Table 1A. Feature dimension summary for the six protocol objects used in this study",
    feature_df,
    header_style_blue,
    title_style_blue
  )
  write_sheet(
    "Network Statistics",
    "Supplemental Table 1B. Network and label statistics for the six protocol objects used in this study",
    network_df,
    header_style_red,
    title_style_red
  )
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)
}

html_escape <- function(text) {
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  gsub(">", "&gt;", text, fixed = TRUE)
}

build_html_table <- function(df, caption, header_color) {
  columns <- names(df)
  header_cells <- paste0("<th>", vapply(columns, html_escape, character(1)), "</th>", collapse = "")
  body_rows <- vapply(seq_len(nrow(df)), function(i) {
    row_class <- if ((i %% 2) == 1) "odd" else "even"
    cells <- paste0(
      ifelse(seq_along(columns) == 1, "<th scope=\"row\">", "<td>"),
      vapply(df[i, ], function(value) html_escape(as.character(value)), character(1)),
      ifelse(seq_along(columns) == 1, "</th>", "</td>"),
      collapse = ""
    )
    sprintf("<tr class=\"%s\">%s</tr>", row_class, cells)
  }, character(1))

  paste0(
    "<div class=\"table-block\">",
    "<table class=\"pub-table\">",
    "<caption>", html_escape(caption), "</caption>",
    "<thead style=\"background:", header_color, ";\"><tr>", header_cells, "</tr></thead>",
    "<tbody>", paste(body_rows, collapse = ""), "</tbody>",
    "</table>",
    "</div>"
  )
}

write_html <- function(feature_df, network_df, output_path) {
  html <- paste0(
    "<!DOCTYPE html>",
    "<html lang=\"en\">",
    "<head>",
    "<meta charset=\"UTF-8\">",
    "<title>Supplemental Table 1 species feature statistics</title>",
    "<style>",
    "body{font-family:Arial,Helvetica,sans-serif;margin:32px;color:#1f2933;background:#ffffff;}",
    "h1{font-size:24px;margin:0 0 8px 0;}",
    "p.lead{margin:0 0 24px 0;color:#52606d;}",
    "section{margin-bottom:36px;}",
    "table.pub-table{border-collapse:collapse;width:100%;max-width:1200px;font-size:14px;}",
    "table.pub-table caption{caption-side:top;text-align:left;font-weight:700;font-size:18px;margin-bottom:10px;color:#102a43;}",
    "table.pub-table th,table.pub-table td{padding:10px 12px;border:1px solid #d9e2ec;text-align:center;}",
    "table.pub-table thead th{color:#ffffff;font-weight:700;}",
    "table.pub-table tbody th{font-weight:700;text-align:left;background:#f8fafc;}",
    "table.pub-table tbody tr.odd td,table.pub-table tbody tr.odd th{background:#ffffff;}",
    "table.pub-table tbody tr.even td,table.pub-table tbody tr.even th{background:#f5f7fa;}",
    ".note{max-width:1200px;font-size:12px;color:#52606d;margin-top:10px;}",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Supplemental Table 1. Species feature statistics</h1>",
    "<p class=\"lead\">All values were derived programmatically from the current frozen-protocol processed feature files and label/split manifests for the six protocol objects used in this study.</p>",
    "<section>",
    build_html_table(
      feature_df,
      "Supplemental Table 1A. Feature dimension summary for the six protocol objects used in this study",
      "#4E79A7"
    ),
    "</section>",
    "<section>",
    build_html_table(
      as.data.frame(lapply(network_df, function(col) if (is.numeric(col)) format_count(col) else col), check.names = FALSE),
      "Supplemental Table 1B. Network and label statistics for the six protocol objects used in this study",
      "#E15759"
    ),
    "<div class=\"note\">Both Fusarium label regimes are reported separately. Feature dimensions were computed from the processed ortholog / expression / localization matrices. Network and label statistics were computed from the processed STRING edge list plus frozen-protocol label and split manifests.</div>",
    "</section>",
    "</body>",
    "</html>"
  )
  writeLines(html, output_path, useBytes = TRUE)
}

table_to_markdown <- function(df) {
  headers <- names(df)
  align <- rep("---", length(headers))
  rows <- vapply(seq_len(nrow(df)), function(i) paste(vapply(df[i, ], as.character, character(1)), collapse = " | "), character(1))
  paste(
    paste(headers, collapse = " | "),
    paste(align, collapse = " | "),
    paste(rows, collapse = "\n"),
    sep = "\n"
  )
}

write_markdown <- function(feature_df, network_df, info_list, output_path) {
  lines <- c(
    "# Supplemental Table 1 species feature statistics",
    "",
    "This file follows the legacy four-species feature-table organization and expands it to the six protocol objects used in the current manuscript.",
    "",
    "## Supplemental Table 1A. Feature dimension summary for the six protocol objects used in this study",
    "",
    table_to_markdown(as.data.frame(feature_df, check.names = FALSE)),
    "",
    "## Supplemental Table 1B. Network and label statistics for the six protocol objects used in this study",
    "",
    table_to_markdown(as.data.frame(network_df, check.names = FALSE)),
    "",
    "## Source audit",
    ""
  )

  for (object_id in object_order) {
    info <- info_list[[object_id]]
    lines <- c(
      lines,
      sprintf("### %s", info$display_name),
      "",
      sprintf("- protocol id: `%s`", object_id),
      sprintf("- orthologs: `%s`", info$files$orthologs),
      sprintf("- expression: `%s`", info$files$expression),
      sprintf("- sublocalization: `%s`", info$files$sublocalization),
      sprintf("- ppi: `%s`", info$files$ppi),
      sprintf("- labels: `%s`", info$files$labels),
      sprintf("- splits: `%s`", info$files$splits),
      ""
    )
  }

  writeLines(lines, output_path, useBytes = TRUE)
}

md_path <- file.path(output_dir, "Supplemental_table1_species_feature.md")
xlsx_path <- file.path(output_dir, "Supplemental_table1_species_feature.xlsx")
html_path <- file.path(output_dir, "Supplemental_table1_species_feature.html")

write_excel(feature_table, network_table, xlsx_path)
write_html(feature_table, network_table, html_path)
write_markdown(feature_table, network_table, dataset_info, md_path)

message("Wrote Markdown: ", md_path)
message("Wrote Excel: ", xlsx_path)
message("Wrote HTML: ", html_path)
