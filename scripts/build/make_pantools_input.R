suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## prepare input data for pangenome construction using PanTools

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

panConf <- confs$data$pangenomes$pecto_demo

cols_metadata <- c(
  "sampleName", "AssemblyAccession",	"AssemblyName", "SpeciesName", "taxonomy_check_status",
  "BioSampleAccn", "BioprojectAccn",
  "strain", "virulence", "virulence_pcr", "geo_loc_country", "continent",
  "host", "isolation_source", "collection_year", "collected_by", "env_broad_scale",
  "type_material", "virulence", "virulence_pcr",
  "source", "type_material", "representative_status", "sample_type",
  "length", "N50", "L50", "n_contigs")

#####################################################################
! dir.exists(panConf$dir) && dir.create(path = panConf$dir, recursive = TRUE)

if(! dir.exists(panConf$dir)){
  dir.create(path = panConf$dir, recursive = TRUE)
}

metadata <- suppressMessages(
  readr::read_tsv(file = confs$data$reference_data$files$metadata)
)

samples <- suppressMessages(
  readr::read_tsv(file = panConf$files$samples)
) %>%
  dplyr::left_join(y = metadata, by = "sampleId") %>%
  dplyr::mutate(
    Genome = 1:n(),
    genomeId = paste("g_", Genome, sep = ""),
    fasta = paste(
      confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""
    ),
    gff3 = paste(
      confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".gff3", sep = ""
    ),
    interpro = paste(
      confs$data$interproscan$dir, "/", sampleId, ".interProScan.gff3", sep = ""
    ),
    cog = paste(
      confs$data$cog$dir, "/", sampleId, ".emapper.annotations", sep = ""
    ),
    ## replace "," with ";"
    dplyr::across(
      .cols = everything(),
      .fns = ~stringr::str_replace_all(
        string = .x, pattern = "[,;]", replacement = " and"
      )
    )
  )

#####################################################################
## FASTA file paths
dplyr::select(samples, fasta) %>%
  readr::write_tsv(
    file = panConf$files$genomes,
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(samples, Genome, gff3) %>%
  readr::write_tsv(file = panConf$files$gff, col_names = FALSE)

## functional annotation file paths
dplyr::select(samples, Genome, interpro) %>%
  readr::write_delim(file = panConf$files$annotations, col_names = FALSE)

## COG annotation file path
dplyr::select(samples, Genome, cog) %>%
  readr::write_delim(file = panConf$files$cog, col_names = FALSE)

## metadata file
dplyr::select(samples, Genome, genomeId, sampleId, !!!cols_metadata) %>%
  readr::write_csv(file = panConf$files$metadata, col_names = TRUE)

## input genome lock file: this file should never change
readr::write_tsv(
  x = dplyr::select(samples, Genome, sampleId),
  file = panConf$files$input_lock
)

## write metadata to excel
wb <- openxlsx::createWorkbook()
currentSheet <- "metadata"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, withFilter = TRUE, keepNA = TRUE, na.string = "NA",
  x = dplyr::select(
    samples, Genome=genomeId, sampleId, !!!cols_metadata
  )
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = file.path(panConf$dir, "pangenome_metadata.xlsx"), overwrite = TRUE
)

