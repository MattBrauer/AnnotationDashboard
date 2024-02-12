library(tidyverse)
library(qs)
library(mazer)
library(ckanr)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(TissueEnrich)
library(shiny)
library(shinydashboard)

base_dir <- here::here()
data_dir <- fs::path(base_dir, "data")
output_dir <- fs::path(base_dir, "output")

ckan_url <- "data.mazetx.com"
ckan_key <- "2975099d-2f83-4935-880d-6fcc8db8b308" # change to environment var
ckanr::ckanr_setup(url = ckan_url, key = ckan_key)
data_science_org_id <- "c121e051-ef65-42d2-a491-62ee1c8a9e9d"
package_id <- '3131be29-48c8-4662-8e3a-6ed4fde7c276'

gene_lists <- qs::qread(fs::path(output_dir, "gene_list.qs"))

gs <- GeneSet(geneIds = gene_lists %>% dplyr::pull(symbol),
              organism = "Homo Sapiens",
              geneIdType = SymbolIdentifier())

output <- teEnrichment(inputGenes = gs)

genePatterns_df <- lapply(names(output[[3]]), function(tissue) {
  df <- assay(output[[3]][[tissue]]) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(tissue = tissue)
  if(nrow(df) == 0 ) df <- tibble::tibble_row(Gene = NA, Group = NA, tissue = tissue)
  return(df)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::bind_rows(tibble::tibble(Gene = geneIds(output[[4]]),
                                  Group = "Unmapped", tissue = NA)) %>%
  dplyr::rename(gene = Gene, group = Group)

genePatterns_list <- setNames(genePatterns_df %>%
                                dplyr::mutate(group = ifelse(stringr::str_ends(group, "Enriched"),
                                                             "enriched",
                                                             "enhanced")) %>%
                                dplyr::mutate(group = ifelse(is.na(tissue), "unmapped", group)) %>%
                                dplyr::group_by(group) %>%
                                dplyr::group_split(),
                              genePatterns_df %>%
                                dplyr::mutate(group = ifelse(stringr::str_ends(group, "Enriched"),
                                                             "enriched",
                                                             "enhanced")) %>%
                                dplyr::mutate(group = ifelse(is.na(tissue), "unmapped", group)) %>%
                                dplyr::group_by(group) %>%
                                dplyr::group_keys() %>%
                                dplyr::pull(group))


genePatterns_comb_mat <- lapply(genePatterns_list,
                                function(class) {
                                  combinations <- class %>%
                                    dplyr::mutate(tissue = ifelse(is.na(tissue), group, tissue),
                                                  membership = 1) %>%
                                    tidyr::pivot_wider(id_cols = gene,
                                                       names_from = tissue,
                                                       values_from = membership,
                                                       values_fill = 0)
                                  comb_mat <- ComplexHeatmap::make_comb_mat(combinations)
                                  return(list(comb = combinations,
                                              comb_mat = comb_mat))
                                })

ht1 = Heatmap(genePatterns_comb_mat[["enriched"]]$comb_mat, name = "Enriched",
              show_row_names = TRUE, show_column_names = TRUE)
ht1 = draw(ht1)

#htShiny(ht1)

up1 = UpSet(genePatterns_comb_mat[["enriched"]]$comb_mat, name = "Enriched",
              show_row_names = TRUE, show_column_names = TRUE)
up1 = draw(up1)

ec <- "enriched"
#ec <- "enhanced"

ui = fluidPage(
  titlePanel("Tissue Expression"),
  fluidRow(

    column(6,
           plotOutput("heatmap", width = 1000, height = 1000, click = "heatmap_click")),
    column(3,
           box(title = uiOutput('title1'),verbatimTextOutput("table1"))),
    column(3,
           box(title = uiOutput('title2'),verbatimTextOutput("table2")))
  )
)

server = function(input, output, session) {

  ht_obj = reactiveVal(NULL)
  ht_pos_obj = reactiveVal(NULL)


  output$heatmap = renderPlot({
    ht1 = draw(UpSet(genePatterns_comb_mat[[ec]]$comb_mat, name = "Enriched",
                     show_row_names = TRUE, show_column_names = TRUE))
    ht1_pos = htPositionsOnDevice(ht1)

    ht_obj(ht1)
    ht_pos_obj(ht1_pos)

  })

  observeEvent(input$heatmap_click, {
    pos = getPositionFromClick(input$heatmap_click)

    selection = selectPosition(ht_obj(), pos, mark = FALSE, ht_pos = ht_pos_obj(),
                               verbose = FALSE)

      row <- selection$row_index[1]
      col <- selection$column_index[1]
      tissue_label <- selection$row_label[1]

      combination_name <- comb_name(genePatterns_comb_mat[[ec]]$comb_mat[col])
      combination_genes <- extract_comb(genePatterns_comb_mat[[ec]]$comb_mat,
                                        combination_name)


      row_genes <- genePatterns_list[[ec]] %>%
        dplyr::filter(tissue == tissue_label) %>%
        dplyr::select(gene) %>%
        dplyr::distinct() %>%
        dplyr::arrange(gene)

      col_genes <- genePatterns_comb_mat[[ec]]$comb[combination_genes,] %>%
        dplyr::select(gene) %>%
        dplyr::left_join(genePatterns_df %>%
                           dplyr::filter(stringr::str_ends(group, "-Enriched")),
                         by="gene") %>%
        dplyr::select(-group) %>%
        dplyr::arrange(gene)

      col_tissues <- col_genes %>%
        dplyr::select(tissue) %>%
        dplyr::distinct() %>%
        dplyr::arrange(tissue) %>%
        dplyr::pull(tissue)

      col_genes <- col_genes %>%
        dplyr::select(gene) %>%
        dplyr::arrange(gene) %>%
        dplyr::distinct()

      output$table1 <- renderPrint(
        cat(row_genes %>% dplyr::pull(gene), sep = "\n")
      )
      output$table2 <- renderPrint(
        cat(col_genes %>% dplyr::pull(gene), sep = "\n")
      )

      output$title1 <- renderUI({
        paste("Union:", tissue_label)
      })

      output$title2 <- renderUI({
        paste("Intersection:\n", paste(col_tissues, collapse = ","))
      })

  })
}
shinyApp(ui, server)
