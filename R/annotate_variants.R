#' Annotate Variants
#'
#' Wrapper function for oncotate from maftools package to facilitate compatibility with VRanges.
#' Returns a dataframe with mutations annotated by Oncotator - requires internet connection. Only compatible with hg19 genome.
#'
#' @param variant_calls \code{VRanges} object from call_all_variants output
#' @importFrom dplyr "%>%"
#' @export
#' @examples
#' \dontrun{
#' annotated_variants <- variant_calls %>% annotate_variants()
#' }
#' @return This function returns a \code{DataFrame} object with annotations for each mutation

annotate_variants <-
function(variant_calls) {

  # Check that input is VRanges
  if(class(variant_calls) != "VRanges"){ stop('Input "variant_calls" needs to be VRanges object') }

  # Check that genome is hg19
  if(!(GenomeInfoDb::genome(variant_calls) %in% c("hg19", "Hg19", "HG19", "GRCh37", "GRCH37", "grch37"))){ stop('Genome needs to be hg19') }

  # get annotation of unique alleles in variant_calls
  anno_full <- variant_calls %>%
    dplyr::as_tibble() %>%
    dplyr::select("chr" = seqnames, start, end, "ref_allele" = ref, "alt_allele" = alt) %>%
    unique() %>%
    maftools::oncotate()

  # keep the relevant columns only
  anno_slim <- anno_full %>%
    dplyr::mutate(chr = paste0("chr", chr),
                  start = as.integer(start),
                  end = as.integer(end)) %>%
    dplyr::select(chr, start, end, "ref" = ref_allele, "alt" = alt_allele,
                  Hugo_Symbol, Strand, Variant_Classification, dbSNP_RS,
                  transcript_id, transcript_strand, transcript_exon, transcript_change, codon_change, protein_change,
                  HGVS_genomic_change, HGVS_coding_DNA_change, dbNSFP_aapos_SIFT, dbNSFP_SIFT_pred, dbNSFP_SIFT_score,
                  `HGNC_Approved Name`, HGNC_Chromosome, `HGNC_Gene family description`, COSMIC_Tissue_tissue_types_affected,
                  COSMIC_Tissue_total_alterations_in_gene, COSMIC_n_overlapping_mutations, COSMIC_overlapping_mutation_AAs,
                  COSMIC_overlapping_primary_sites, `CGC_Cancer Germline Mut`, `CGC_Cancer Somatic Mut`,
                  `CGC_Cancer Syndrome`, `CGC_Tumour Types  (Somatic Mutations)`, `CGC_Tumour Types (Germline Mutations)`,
                  UniProt_GO_Biological_Process, UniProt_GO_Cellular_Component, UniProt_GO_Molecular_Function)

  # left join samp with annotation
  annotated <- variant_calls %>%
    dplyr::as_tibble() %>%
    dplyr::rename("chr" = seqnames) %>%
    dplyr::mutate(chr = as.character(chr)) %>%
    dplyr::left_join(., anno_slim, by = c("chr", "start", "end", "ref", "alt"))

  return(annotated)
}
