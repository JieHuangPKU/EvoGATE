#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(readr)
  library(yaml)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(project_root, "results", "tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_order <- c("celegans", "scerevisiae", "fgraminearum", "human")
species_labels <- c(
  celegans = "C. elegans",
  scerevisiae = "S. cerevisiae",
  fgraminearum = "F. graminearum",
  human = "H. sapiens"
)
feature_rows <- c(
  orthologs = "Ortholog",
  expression = "Expression",
  sublocalization = "Localization"
)
config_paths <- c(
  human = "configs/epgat_graph_benchmark_human.yaml",
  celegans = "configs/epgat_graph_benchmark_celegans.yaml",
  scerevisiae = "configs/epgat_graph_benchmark_scerevisiae.yaml",
  fgraminearum = "configs/epgat_graph_benchmark_fgraminearum.yaml"
)

format_count <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

run_python <- function(code, args = character()) {
  script_path <- tempfile(fileext = ".py")
  writeLines(code, script_path, useBytes = TRUE)
  on.exit(unlink(script_path), add = TRUE)
  output <- system2("python", c(script_path, args), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) {
    stop(paste(output, collapse = "\n"))
  }
  output
}

read_npy_shape <- function(path) {
  code <- paste(
    "import numpy as np, sys",
    "arr = np.load(sys.argv[1], allow_pickle=True)",
    "print('\\t'.join(str(x) for x in arr.shape))",
    sep = "\n"
  )
  shape_text <- run_python(code, normalizePath(path, winslash = "/", mustWork = TRUE))
  as.integer(strsplit(shape_text[[1]], "\t", fixed = TRUE)[[1]])
}

read_dataset_info <- function(species_key, config_path) {
  cfg <- yaml::read_yaml(file.path(project_root, config_path))
  dataset_dir <- file.path(project_root, cfg$dataset$source_dir)

  feature_schema_path <- file.path(dataset_dir, "feature_schema.tsv")
  feature_matrix_path <- file.path(dataset_dir, "feature_matrix.npy")
  node_manifest_path <- file.path(dataset_dir, "node_manifest.tsv")
  edge_index_path <- file.path(dataset_dir, "edge_index.npy")
  edge_table_path <- file.path(dataset_dir, "edge_table.tsv")
  label_manifest_path <- file.path(dataset_dir, "label_manifest.tsv")

  feature_schema <- readr::read_tsv(feature_schema_path, show_col_types = FALSE)
  node_manifest <- readr::read_tsv(node_manifest_path, show_col_types = FALSE)
  edge_table <- readr::read_tsv(edge_table_path, show_col_types = FALSE)
  label_manifest <- readr::read_tsv(
    label_manifest_path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x), "", .x)))

  feature_matrix_shape <- read_npy_shape(feature_matrix_path)
  edge_index_shape <- read_npy_shape(edge_index_path)

  if (length(feature_matrix_shape) != 2) {
    stop(sprintf("Unexpected feature_matrix shape for %s: %s", species_key, paste(feature_matrix_shape, collapse = "x")))
  }
  if (length(edge_index_shape) != 2 || edge_index_shape[[2]] != 2) {
    stop(sprintf("Unexpected edge_index shape for %s: %s", species_key, paste(edge_index_shape, collapse = "x")))
  }

  schema_dimension_total <- sum(feature_schema$dimension)
  if (schema_dimension_total != feature_matrix_shape[[2]]) {
    stop(
      sprintf(
        "Feature schema total (%s) does not match feature matrix columns (%s) for %s",
        schema_dimension_total,
        feature_matrix_shape[[2]],
        species_key
      )
    )
  }
  if (nrow(node_manifest) != feature_matrix_shape[[1]]) {
    stop(
      sprintf(
        "Node manifest rows (%s) do not match feature matrix rows (%s) for %s",
        nrow(node_manifest),
        feature_matrix_shape[[1]],
        species_key
      )
    )
  }
  if (nrow(edge_table) != edge_index_shape[[1]]) {
    stop(
      sprintf(
        "Edge table rows (%s) do not match edge index rows (%s) for %s",
        nrow(edge_table),
        edge_index_shape[[1]],
        species_key
      )
    )
  }

  labeled_manifest <- label_manifest %>%
    mutate(
      is_labeled_flag = tolower(is_labeled) %in% c("true", "1", "yes"),
      label_numeric = suppressWarnings(as.integer(label))
    ) %>%
    filter(is_labeled_flag)

  feature_dimensions <- vapply(names(feature_rows), function(block) {
    matched <- feature_schema %>% filter(feature_block == block)
    if (nrow(matched) == 0) {
      return(0L)
    }
    as.integer(sum(matched$dimension))
  }, integer(1))

  list(
    species = species_key,
    species_label = unname(species_labels[[species_key]]),
    config_path = config_path,
    dataset_dir = cfg$dataset$source_dir,
    files = list(
      feature_schema = sub(paste0("^", project_root, "/?"), "", feature_schema_path),
      feature_matrix = sub(paste0("^", project_root, "/?"), "", feature_matrix_path),
      node_manifest = sub(paste0("^", project_root, "/?"), "", node_manifest_path),
      edge_index = sub(paste0("^", project_root, "/?"), "", edge_index_path),
      edge_table = sub(paste0("^", project_root, "/?"), "", edge_table_path),
      label_manifest = sub(paste0("^", project_root, "/?"), "", label_manifest_path)
    ),
    feature_dimensions = feature_dimensions,
    node_count = nrow(node_manifest),
    edge_count = edge_index_shape[[1]],
    labeled_count = nrow(labeled_manifest),
    essential_count = sum(labeled_manifest$label_numeric == 1, na.rm = TRUE),
    test_count = sum(labeled_manifest$split == "test", na.rm = TRUE)
  )
}

dataset_info <- lapply(species_order, function(species_key) {
  read_dataset_info(species_key, config_paths[[species_key]])
})
names(dataset_info) <- species_order

feature_table <- data.frame("Feature Type" = unname(feature_rows), check.names = FALSE)
for (species_key in species_order) {
  feature_table[[species_labels[[species_key]]]] <- unname(dataset_info[[species_key]]$feature_dimensions[names(feature_rows)])
}

network_table <- data.frame(
  Statistic = c("N.nodes", "N.edges", "N.labeled genes", "N.essential genes", "N.test genes"),
  check.names = FALSE
)
for (species_key in species_order) {
  info <- dataset_info[[species_key]]
  network_table[[species_labels[[species_key]]]] <- c(
    info$node_count,
    info$edge_count,
    info$labeled_count,
    info$essential_count,
    info$test_count
  )
}

write_excel <- function(feature_df, network_df, output_path) {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Feature Dimensions", gridLines = FALSE)
  openxlsx::addWorksheet(wb, "Network Statistics", gridLines = FALSE)

  title_style_blue <- openxlsx::createStyle(
    fontSize = 14, textDecoration = "bold", halign = "left"
  )
  title_style_red <- openxlsx::createStyle(
    fontSize = 14, textDecoration = "bold", halign = "left"
  )
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
      openxlsx::addStyle(
        wb, sheet_name, fill_style,
        rows = row_index, cols = seq_len(ncol(df)),
        gridExpand = TRUE, stack = TRUE
      )
    }
    openxlsx::setColWidths(wb, sheet_name, cols = seq_len(ncol(df)), widths = "auto")
    openxlsx::freezePane(wb, sheet_name, firstActiveRow = 4)
  }

  write_sheet(
    "Feature Dimensions",
    "Table 1. Feature dimension summary for the four species used in this study",
    feature_df,
    header_style_blue,
    title_style_blue
  )
  write_sheet(
    "Network Statistics",
    "Table 2. Network and label statistics for the four species used in this study",
    network_df,
    header_style_red,
    title_style_red
  )
  openxlsx::saveWorkbook(wb, output_path, overwrite = TRUE)

  sanitize_code <- paste(
    "import os",
    "import sys",
    "import tempfile",
    "import zipfile",
    "from xml.etree import ElementTree as ET",
    "",
    "src = sys.argv[1]",
    "fd, tmp = tempfile.mkstemp(suffix='.xlsx', dir=os.path.dirname(src) or '.')",
    "os.close(fd)",
    "",
    "with zipfile.ZipFile(src, 'r') as zin, zipfile.ZipFile(tmp, 'w', zipfile.ZIP_DEFLATED) as zout:",
    "    for item in zin.infolist():",
    "        data = zin.read(item.filename)",
    "        if item.filename.startswith('xl/worksheets/_rels/') and item.filename.endswith('.rels'):",
    "            root = ET.fromstring(data)",
    "            kept = []",
    "            for rel in root.findall('{http://schemas.openxmlformats.org/package/2006/relationships}Relationship'):",
    "                rel_type = rel.attrib.get('Type', '')",
    "                if rel_type.endswith('/drawing') or rel_type.endswith('/vmlDrawing'):",
    "                    continue",
    "                kept.append(rel)",
    "            if kept:",
    "                new_root = ET.Element(root.tag, root.attrib)",
    "                for rel in kept:",
    "                    new_root.append(rel)",
    "                data = ET.tostring(new_root, encoding='utf-8', xml_declaration=True)",
    "            else:",
    "                continue",
    "        zout.writestr(item, data)",
    "",
    "os.replace(tmp, src)",
    sep = "\n"
  )
  invisible(run_python(sanitize_code, normalizePath(output_path, winslash = "/", mustWork = TRUE)))
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
    "<title>Four-species feature statistics</title>",
    "<style>",
    "body{font-family:Arial,Helvetica,sans-serif;margin:32px;color:#1f2933;background:#ffffff;}",
    "h1{font-size:24px;margin:0 0 8px 0;}",
    "p.lead{margin:0 0 24px 0;color:#52606d;}",
    "section{margin-bottom:36px;}",
    "table.pub-table{border-collapse:collapse;width:100%;max-width:960px;font-size:14px;}",
    "table.pub-table caption{caption-side:top;text-align:left;font-weight:700;font-size:18px;margin-bottom:10px;color:#102a43;}",
    "table.pub-table th,table.pub-table td{padding:10px 12px;border:1px solid #d9e2ec;text-align:center;}",
    "table.pub-table thead th{color:#ffffff;font-weight:700;}",
    "table.pub-table tbody th{font-weight:700;text-align:left;background:#f8fafc;}",
    "table.pub-table tbody tr.odd td,table.pub-table tbody tr.odd th{background:#ffffff;}",
    "table.pub-table tbody tr.even td,table.pub-table tbody tr.even th{background:#f5f7fa;}",
    ".note{max-width:960px;font-size:12px;color:#52606d;margin-top:10px;}",
    "</style>",
    "</head>",
    "<body>",
    "<h1>Four-species feature statistics</h1>",
    "<p class=\"lead\">All values were derived programmatically from the current benchmark dataset directories used by the archived four-species graph benchmark pipeline.</p>",
    "<section>",
    build_html_table(
      feature_df,
      "Table 1. Feature dimension summary for the four species used in this study",
      "#4E79A7"
    ),
    "</section>",
    "<section>",
    build_html_table(
      as.data.frame(lapply(network_df, function(col) if (is.numeric(col)) format_count(col) else col), check.names = FALSE),
      "Table 2. Network and label statistics for the four species used in this study",
      "#E15759"
    ),
    "<div class=\"note\">Feature dimensions were taken from <code>feature_schema.tsv</code> and validated against <code>feature_matrix.npy</code>. Network and label statistics were derived from <code>node_manifest.tsv</code>, <code>edge_index.npy</code>, and <code>label_manifest.tsv</code>, with <code>edge_table.tsv</code> used as a row-count cross-check for the graph edge list.</div>",
    "</section>",
    "</body>",
    "</html>"
  )
  writeLines(html, output_path, useBytes = TRUE)
}

write_audit_readme <- function(info_list, output_path) {
  lines <- c(
    "# Four-species feature table audit note",
    "",
    "This note documents the authoritative inputs used to build the publication-style four-species descriptive tables.",
    "",
    "## Scope",
    "",
    "The current archived four-species benchmark runner is `src/train/run_epgat_graph_benchmark.py`. It reads these benchmark configs:",
    ""
  )

  for (species_key in species_order) {
    info <- info_list[[species_key]]
    lines <- c(
      lines,
      sprintf("- `%s` -> `%s` -> dataset directory `%s`", species_key, info$config_path, info$dataset_dir)
    )
  }

  lines <- c(
    lines,
    "",
    "The dataset directories above were treated as authoritative because they are the exact `dataset.source_dir` values consumed by `src/train/train_epgat_graph_models.py` for the current four-species replay benchmark.",
    "",
    "## Authoritative files and why they were chosen",
    "",
    "- `feature_schema.tsv`: authoritative source for feature block membership and per-block dimensions because the trainer copies this schema into each model run and it explicitly records the `feature_block` and `dimension` columns used to describe the model input layout.",
    "- `feature_matrix.npy`: authoritative validation source for the final feature matrix shape because it is the actual matrix loaded by the trainer. Its column count was checked against the summed `feature_schema.tsv` dimensions and its row count was checked against `node_manifest.tsv`.",
    "- `node_manifest.tsv`: authoritative source for graph node count because it is the node universe aligned to the model input and is read directly by the trainer.",
    "- `edge_index.npy`: authoritative source for graph edge count because it is the graph edge object loaded directly by the trainer for model execution.",
    "- `edge_table.tsv`: cross-check for `edge_index.npy` row count because it is the tabular representation of the same graph edges in the dataset directory.",
    "- `label_manifest.tsv`: authoritative source for label availability, positive labels, and split assignments because the trainer derives `labeled`, `train`, `val`, and `test` subsets directly from this file.",
    "",
    "## How each statistic was computed",
    "",
    "### Table 1",
    "",
    "- `Ortholog`: sum of `dimension` values in `feature_schema.tsv` where `feature_block == \"orthologs\"`.",
    "- `Expression`: sum of `dimension` values in `feature_schema.tsv` where `feature_block == \"expression\"`.",
    "- `Localization`: sum of `dimension` values in `feature_schema.tsv` where `feature_block == \"sublocalization\"`.",
    "- Validation: the total dimension across all rows in `feature_schema.tsv` must equal the number of columns in `feature_matrix.npy`.",
    "",
    "### Table 2",
    "",
    "- `N.nodes`: number of rows in `node_manifest.tsv`. This was validated against the number of rows in `feature_matrix.npy`.",
    "- `N.edges`: number of rows in `edge_index.npy`. This was validated against the number of rows in `edge_table.tsv`.",
    "- `N.labeled genes`: number of rows in `label_manifest.tsv` where `is_labeled` is one of `True`, `true`, `1`, or `yes`.",
    "- `N.essential genes`: number of labeled rows in `label_manifest.tsv` where `label == 1`.",
    "- `N.test genes`: number of labeled rows in `label_manifest.tsv` where `split == \"test\"`.",
    "",
    "## Species-specific dataset files",
    ""
  )

  for (species_key in species_order) {
    info <- info_list[[species_key]]
    lines <- c(
      lines,
      sprintf("### %s", info$species_label),
      "",
      sprintf("- `feature_schema.tsv`: `%s`", info$files$feature_schema),
      sprintf("- `feature_matrix.npy`: `%s`", info$files$feature_matrix),
      sprintf("- `node_manifest.tsv`: `%s`", info$files$node_manifest),
      sprintf("- `edge_index.npy`: `%s`", info$files$edge_index),
      sprintf("- `edge_table.tsv`: `%s`", info$files$edge_table),
      sprintf("- `label_manifest.tsv`: `%s`", info$files$label_manifest),
      ""
    )
  }

  writeLines(lines, output_path, useBytes = TRUE)
}

excel_path <- file.path(output_dir, "four_species_feature_statistics.xlsx")
html_path <- file.path(output_dir, "four_species_feature_statistics.html")
readme_path <- file.path(output_dir, "README_four_species_feature_tables.md")

write_excel(feature_table, network_table, excel_path)
write_html(feature_table, network_table, html_path)
write_audit_readme(dataset_info, readme_path)

message("Wrote Excel workbook: ", excel_path)
message("Wrote HTML report: ", html_path)
message("Wrote audit note: ", readme_path)
