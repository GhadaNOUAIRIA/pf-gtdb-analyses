#!/usr/bin/env Rscript

# select_sequences_in_GTDB_representatives.R
#
# Subsets sequences from genomes in the GTDB representatives list.
# Writes to feather files.
#
# Author: daniel.lundin@dbb.su.se, nouairia.ghada@gmail.com
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

accessions <- read_feather('pfitmap-gtdb.accessions.feather')
domains <- read_feather('pfitmap-gtdb.domains.feather')
proteins <- read_feather('pfitmap-gtdb.proteins.feather')
sequences <- read_feather('pfitmap-gtdb.sequences.feather')
taxa <- read_feather('pfitmap-gtdb.taxa.feather')
tblout <- read_feather('pfitmap-gtdb.tblout.feather')
domtblout <- read_feather('pfitmap-gtdb.domtblout.feather')

representatives <- read_tsv('sp_clusters.tsv') 
names(representatives)<-str_replace_all(names(representatives), c(" " = "."))

Representatives_taxa = taxa %>%
  mutate(genome = sub("GCA_", "", genome_accno)) %>%
  semi_join(representatives %>% 
              mutate (genome = sub("GB_GCA_", "", Representative.genome)) %>% mutate(genome = sub("RS_GCF_", "", genome)) %>% 
              mutate(genome = sub("\\.[1-9]", "", genome)), 
            by = 'genome')

Representatives_taxa %>%
  select(genome_accno, tdomain, tphylum, tclass, torder, tfamily, tgenus, tspecies, trank, ncbi_taxon_id) %>%
  write_feather("pfitmap-gtdb-rep.taxa.feather")

Representatives_accessions = accessions %>%
  inner_join(Representatives_taxa, by = 'genome_accno') %>%
  select(accno)

accessions %>%
  inner_join (Representatives_accessions, by = 'accno') %>%
  write_feather("pfitmap-gtdb-rep.accessions.feather")

proteins %>% 
    inner_join (Representatives_accessions, by = 'accno') %>%
    write_feather("pfitmap-gtdb-rep.proteins.feather")

domains %>% 
  inner_join (Representatives_accessions, by = 'accno') %>%
  write_feather("pfitmap-gtdb-rep.domains.feather")

sequences %>%
  inner_join (Representatives_accessions, by = 'accno') %>%
  write_feather("pfitmap-gtdb-rep.sequences.feather")

tblout %>% 
  inner_join (Representatives_accessions, by = 'accno') %>%
  write_feather("pfitmap-gtdb-rep.tblout.feather")

domtblout %>%
  inner_join (Representatives_accessions, by = 'accno') %>%
  write_feather("pfitmap-gtdb-rep.domtblout.feather")