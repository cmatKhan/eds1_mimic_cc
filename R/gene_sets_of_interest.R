library(tidyverse)
library(here)

combineSgdManualAndCompGOgenes = function(manual_path, comp_path){

  manual_df = read_tsv(here(manual_path), skip = 8)
  comp_df = read_tsv(here(comp_path), skip = 8)

  bind_rows(manual_df, comp_df) %>%
    distinct(`Gene/Complex`, `Systematic Name/Complex Accession`) %>%
    dplyr::rename(symbol = `Gene/Complex`,
                  id = `Systematic Name/Complex Accession`)
}

set_paths = list(
  lysine_genes = c(here("data/lysine_biosynthetic_process_via_aminoadipic_acid_manual.txt"),
                   here("data/lysine_biosynthetic_process_via_aminoadipic_acid_comp.txt")),
  aerobic_resp_genes = c(here("data/aerobic_respiration_manual.txt"),
                         here("data/aerobic_respiration_comp.txt")),
  anaerobic_resp_genes = c(here("data/anaerobic_respiration_annotations.txt"),
                           here("data/anaerobic_respiration_annotations.txt")),
  tca_genes = c(here("data/tricarboxylic_acid_cycle_manual.txt"),
                here("data/tricarboxylic_acid_cycle_comp.txt")),
  glyoxylate_genes = c(here("data/glyoxylate_cycle_manual.txt"),
                       here("data/glyoxylate_cycle_comp.txt")),
  glycolytic_genes = c(here("data/glycolytic_process_manual.txt"),
                       here("data/glycolytic_process_comp.txt")),
  membrane_transport_genes = c(here("data/transmembrane_transport_manual.txt"),
                               here("data/transmembrane_transport_comp.txt")),
  response_to_glucose_genes = c(here("data/response_to_glucose_annotations.txt"),
                                here("data/response_to_glucose_annotations.txt")),
  hexose_transmem_trans_genes = c(here("data/hexose_transmembrane_transport_manual.txt"),
                                  here("data/hexose_transmembrane_transport_manual.txt")),
  glucose_transmem_trans_genes = c(here("data/glucose_transmembrane_transport_manual.txt"),
                                   here("data/glucose_transmembrane_transport_comp.txt")),
  glucose_transmem_transport_activity_genes = c(here("data/glucose_transmembrane_transporter_activity_annotations.txt"),
                                                here("data/glucose_transmembrane_transporter_activity_annotations.txt")),
  glucose_import_genes = c(here("data/glucose_import_annotations.txt"),
                           here("data/glucose_import_annotations.txt")),
  hexokinase_genes = c(here("data/hexokinase_activity_manual.txt"),
                       here("data/hexokinase_activity_comp.txt")),
  acid_phosphate_activity = c(here("data/acid_phosphatase_activity_manual.txt"),
                              here("data/acid_phosphatase_activity_manual.txt")),
  phoscontaining_compound_meta_proc = c(here("data/phosphatecontaining_compound_metabolic_process_manual.txt"),
                                        here("data/phosphatecontaining_compound_metabolic_process_comp.txt")),
  phosphate_ion_transmem_transport = c(here("data/phosphate_ion_transmembrane_transport_manual.txt"),
<<<<<<< HEAD
                                       here("data/phosphate_ion_transmembrane_transport_comp.txt")),
  cellular_amino_acid_metabolic_process = c(here("data/cellular_amino_acid_metabolic_process_annotations_manual.txt"),
                                       here("data/cellular_amino_acid_metabolic_process_annotations_comp.txt"))
=======
                                       here("data/phosphate_ion_transmembrane_transport_comp.txt"))
>>>>>>> 8fa902069a054d712991ad24f35a0a2007b11ed5

)


gene_set_df_list = map(set_paths, ~combineSgdManualAndCompGOgenes(.[[1]], .[[2]]))
