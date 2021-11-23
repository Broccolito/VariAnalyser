#' Annotate the GWAS results and include pathway and disease information
#'
#' @param variant_file The name of the text file containing GWAS statistics
#' @param chr The name of the chromosome number column in the text file
#' @param pos The name of the chromosome position column in the text file
#' @param ref The name of the reference allele column in the text file
#' @param alt The name of the alternative allele column in the text file
#' @param pvalue The name of the pvalue column in the text file
#' @param delimiter The delimiter used in the text file
#' @param write_excel Whether to write the result to an excel file, else write as CSV
#'
#' @return A data frame of the annotated results
#'
#' @examples
#' PathwayDiseaseAnnotator("data/example/FHS_EA_MRS_merged.txt")


PathwayDiseaseAnnotator = function(variant_file = "data/example/FHS_EA_MRS_5e8_snplist.txt",
                                   chr = "CHR",
                                   pos = "POS",
                                   ref = "NEA",
                                   alt = "EA",
                                   delimiter = " ",
                                   write_excel = TRUE){

  data("genes")
  data("diseases_pathways_flatten")
  data("gene_pathways_flatten")

  cat("Annotating variants...\n")
  variants_annotated = VariantAnnotator(variant_file,
                                        chr,
                                        pos,
                                        ref,
                                        alt,
                                        delimiter,
                                        write_csv = FALSE)

  cat("Annotating Diseases and Pathways...\n")
  data_annotated = variants_annotated %>%
    mutate(GeneSymbol = gene_symbol) %>%
    left_join(select(genes, GeneSymbol, GeneID),
              by = "GeneSymbol") %>%
    left_join(gene_pathways_flatten, "GeneSymbol") %>%
    left_join(diseases_pathways_flatten, by = "GeneSymbol") %>%
    select(-GeneSymbol)

  cat("Saving Files...\n")
  if(write_excel){
    cat("Saving as an excel file...\n")
    file_name = gsub(pattern = ".txt",
                     replacement = "_PathwayDiseaseAnnotated.xlsx",
                     variant_file)
    write_xlsx(data_annotated, path = file_name)
  }else{
    cat("Saving as a csv file...\n")
    file_name = gsub(pattern = ".txt",
                     replacement = "_PathwayDiseaseAnnotated.csv",
                     variant_file)
    write.csv(data_annotated, file = file_name,
              quote = FALSE, row.names = FALSE)
  }

  return(data_annotated)

}










