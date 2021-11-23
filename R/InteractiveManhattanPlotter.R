#' Plot an interactive Manhattan plot and save it as an HTML
#'
#' @param gwas_result_file The name of the text file containing GWAS statistics
#' @param chr The name of the chromosome number column in the text file
#' @param pos The name of the chromosome position column in the text file
#' @param ref The name of the reference allele column in the text file
#' @param alt The name of the alternative allele column in the text file
#' @param pvalue The name of the pvalue column in the text file
#' @param delimiter The delimiter used in the text file
#' @param color1 First color of the Manhattan plot
#' @param color2 Second color of the Manhattan plot
#' @param threshold_linecolor Default is set as black
#'
#' @return A html widget of the Manhattan plot generated
#'
#' @examples
#' InteractiveManhattanPlotter("data/example/FHS_EA_MRS_merged.txt")


InteractiveManhattanPlotter = function(gwas_result_file = "data/example/FHS_EA_MRS_merged.txt",
                                       chr = "CHR",
                                       pos = "POS",
                                       ref = "NEA",
                                       alt = "EA",
                                       pvalue = "p.value",
                                       delimiter = " ",
                                       color1 = "cornflowerblue",
                                       color2 = "blue4",
                                       threshold_linecolor = "black"){

  cat("Loading GWAS Statistics...\n")
  data = read.delim(gwas_result_file, sep = delimiter)
  cat("GWAS Statistics Loaded...\n")

  names(data)[which(names(data)==chr)] = "chr"
  names(data)[which(names(data)==pos)] = "pos"
  names(data)[which(names(data)==ref)] = "ref"
  names(data)[which(names(data)==alt)] = "alt"
  names(data)[which(names(data)==pvalue)] = "pvalue"

  data = data %>%
    select(chr, pos, ref, alt, pvalue) %>%
    mutate(chr = ifelse(chr == "X", "23", chr)) %>%
    mutate(chr = as.numeric(chr)) %>%
    mutate(pos = as.numeric(pos)) %>%
    mutate(pvalue = as.numeric(pvalue))

  data_chromosome_position = data %>%
    group_by(chr) %>%
    summarise(chromosome_length = max(pos)) %>%
    mutate(tot = cumsum(as.numeric(chromosome_length)) - chromosome_length) %>%
    select(-chromosome_length)

  data_cleaned = left_join(data, data_chromosome_position, by = "chr") %>%
    arrange(chr, pos) %>%
    mutate(bp_cum = pos + tot) %>%
    mutate(nlogp = -log(pvalue, base = 10)) %>%
    select(chr, pos, ref, alt, pvalue, nlogp, tot, bp_cum)

  axisdf = data_cleaned %>%
    group_by(chr) %>% summarize(center = (max(bp_cum) + min(bp_cum))/2)

  data_cleaned = data_cleaned %>%
    filter(nlogp >= 5) %>%
    mutate(query = paste0("chr", chr, ":g.", pos, ref, ">", alt))

  if(dim(data_cleaned)[1]<1){
    return("No Significant Signal Detected...\n")
  }

  cat("Generating Variant Annotation...\n")
  variant_client = BioThingsClient("variant")
  annotation = vector()
  for(q in data_cleaned$query){
    try({
      v = btGet(variant_client, q)[[1]]

      gene_symbol = NA
      gene_name = NA
      variant_type = NA
      variant_consequence = NA
      variant_CADD_phred = NA
      variant_phylop_vertebrate = NA
      variant_polyphen = NA
      variant_sift = NA
      variant_rs = NA

      gene_symbol = paste(v[["dbsnp"]][["gene"]][["symbol"]], collapse = "; ")
      gene_name = gsub(pattern = ",", ";", v[["dbsnp"]][["gene"]][["name"]])
      variant_type = gsub(",","; ", paste(v[["cadd"]][["consdetail"]],
                                          collapse = "; "))
      variant_consequence = gsub(",","; ", paste(v[["cadd"]][["consequence"]],
                                                 collapse = "; "))
      variant_CADD_phred = v[["cadd"]][["phred"]]
      variant_phylop_vertebrate = v[["cadd"]][["phylop"]][["vertebrate"]]
      variant_polyphen = v[["cadd"]][["polyphen"]][["val"]]
      variant_sift = v[["cadd"]][["sift"]][["val"]]
      variant_rs = v[["dbsnp"]][["rsid"]]

      gene_symbol = ifelse(is.null(gene_symbol), NA, gene_symbol)
      gene_name = ifelse(is.null(gene_name), NA, gene_name)
      variant_type = ifelse(is.null(variant_type), NA, variant_type)
      variant_consequence = ifelse(is.null(variant_consequence), NA, variant_consequence)

      variant_CADD_phred = ifelse(is.null(variant_CADD_phred), NA, variant_CADD_phred)
      variant_phylop_vertebrate = ifelse(is.null(variant_phylop_vertebrate),
                                         NA, variant_phylop_vertebrate)
      variant_polyphen = ifelse(is.null(variant_polyphen), NA, variant_polyphen)
      variant_sift = ifelse(is.null(variant_sift), NA, variant_sift)
      variant_rs = ifelse(is.null(variant_rs), NA, variant_rs)

      anno = tibble(query = q, gene_symbol, gene_name, variant_type,
                    variant_consequence, variant_CADD_phred,
                    variant_phylop_vertebrate, variant_polyphen,
                    variant_sift, variant_rs)

      annotation = rbind.data.frame(annotation, anno)
      cat(paste0("Annotating ", q," ..\n"))
    }, silent = TRUE)
  }

  data_cleaned = data_cleaned %>%
    left_join(annotation, by = "query")

  cat("Generating Interactive Manhattan Plot of significant Signals...\n")
  suppressWarnings({
    plt = ggplot(data_cleaned, aes(x = bp_cum, y = nlogp)) +
      geom_point(aes(color=as.factor(chr),
                     query = query,
                     variant_rs = variant_rs,
                     gene_symbol = gene_symbol,
                     gene_name = gene_name,
                     variant_type = variant_type,
                     variant_consequence = variant_consequence,
                     variant_CADD_phred = variant_CADD_phred,
                     variant_phylop_vertebrate = variant_phylop_vertebrate,
                     variant_polyphen = variant_polyphen,
                     variant_sift = variant_sift),
                 alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c(color1, color2), 22)) +
      scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
      scale_y_continuous(limits = c(4, max(max(data_cleaned$nlogp)+1,
                                           -log(5e-8,base = 10)+1)),
                         expand = c(0, 0)) +
      geom_hline(yintercept = -log(5e-8,base = 10),
                 linetype = "solid", color = threshold_linecolor) +
      geom_hline(yintercept = -log(1e-5,base = 10),
                 linetype = "dashed", color = threshold_linecolor) +
      xlab("Chromosome") +
      ylab("-log(p-value)") +
      theme_pubr() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)
      )
  })

  plt = ggplotly(plt, tooltip = c("nlogp", "query", "variant_rs",
                                  "gene_symbol", "gene_name", "variant_type",
                                  "variant_consequence", "variant_CADD_phred",
                                  "variant_phylop_vertebrate", "variant_polyphen",
                                  "variant_sift"))

  plot_name = gsub(pattern = ".txt",
                   replacement = "_InteractiveManhattanPlot.html",
                   gwas_result_file)
  saveWidget(plt, plot_name, selfcontained = TRUE, libdir = "lib")

  return(plt)

}

