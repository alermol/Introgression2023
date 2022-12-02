# Function return the table containing either table containing number of SNPs of
# different origin in each window on chromosome or vertical mirror bar plot for
# each chromosome representing number of species-specific SNP in each window

count_SNP_by_window <- function(geno_data,
                                window_size = 10e6,
                                chr_data = 'source_data/chr_length.tsv',
                                plot = T) {
  geno_data <- readr::read_tsv(geno_data)
  
  cols <- colnames(geno_data)[4:6]
  
  geno_data <- geno_data |> 
    dplyr::mutate(orig = dplyr::case_when(
      rlang::parse_expr(paste('`',  cols[1], '`', '==', '`', cols[3], '`', 
                              '&', '`', cols[2], '`', '!=', '`', cols[3], 
                              '`', '~ "TA"', sep = '')),
      rlang::parse_expr(paste('`',  cols[1], '`', '!=', '`', cols[3], '`', 
                              '&', '`', cols[2], '`', '==', '`', cols[3], 
                              '`', '~ "TD"', sep = '')),
      T ~ 'UN'
    ))
  
  chr_lengths_data <- read.delim(chr_data, header = F, sep = '\t') |>
    dplyr::mutate(Start = 0,
                  .after = V1) |>
    dplyr::mutate(data = purrr::map(V2, ~ tibble::tibble(
      chunk = seq(0, ..1 / window_size) + 1,
      start = (chunk - 1) * window_size,
      end = start + window_size
    ))) |>
    tidyr::unnest(data) |> 
    dplyr::rename('chr' = 'V1', 'End' = 'V2') |> 
    dplyr::mutate(
      data = purrr::pmap(
        list(chr, start, end),
        ~ {
          geno_data |> 
            dplyr::filter(Chr_id == ..1,
                          dplyr::between(Start, ..2, ..3)) |> 
            dplyr::summarise(
              snp_num_total = nrow(dplyr::cur_data()),
              snp_num_TA = sum(dplyr::cur_data()$orig == 'TA'),
              snp_num_TD = sum(dplyr::cur_data()$orig == 'TD'),
              snp_num_UN = sum(dplyr::cur_data()$orig == 'UN')
          )
        }
      )
    ) |> 
    tidyr::unnest(data)

  if (plot == T) {
    annot_pos <- -max(chr_lengths_data$snp_num_TA) - 2
    
    plot <- ggplot2::ggplot(data = chr_lengths_data,
                    mapping =  ggplot2::aes(x = start,
                                            group = 0)) +
      ggplot2::geom_col(mapping = ggplot2::aes(y = -snp_num_TA, 
                                               fill = ggsci::pal_igv()(2)[1]),
                        color = 'white', position = 'dodge', size = 0.2) +
      ggplot2::geom_col(mapping = ggplot2::aes(y = snp_num_TD,
                                               fill = ggsci::pal_igv()(2)[2]),
                        color = 'white', position = 'dodge', size = 0.2) +
      ggplot2::geom_rect(mapping = ggplot2::aes(
        xmin = -(window_size / 2), 
        ymin = -0.5, 
        xmax = end - window_size / 2, 
        ymax = 0.5
      )) +
      ggplot2::facet_wrap( ~ chr, nrow = 2, ncol = 7, 
                           scales = 'fixed', dir = 'v') +
      ggplot2::labs(title = glue::glue('{cols[1]} x {cols[2]} >  {cols[3]}'),
                    fill = 'Origin of SNPs') +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.text = ggtext::element_markdown(),
                     legend.position = 'right',
                     axis.title = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_line(),
                     strip.text.x = ggplot2::element_text(size = 14),
                     axis.text.y = ggplot2::element_text(size = 12)) +
      ggplot2::scale_fill_manual(name = 'Origin of SNPs',
                                 values = c(ggsci::pal_igv()(2)[1],
                                            ggsci::pal_igv()(2)[2]),
                                 labels = c('*Triticum aestivum*',
                                            '*Triticum durum*')) + 
      ggplot2::coord_flip() + 
      ggplot2::scale_x_reverse(labels = function(x) glue::glue('{round(x / 1e6, 0)} Mbp'),
                               breaks = function(x) seq(0, max(x), by = 100e6)) + 
      ggplot2::scale_y_continuous(labels = function(x) stringr::str_remove(x, '-'),
                                  limits = c(-max(c(chr_lengths_data$snp_num_TA, 
                                                    chr_lengths_data$snp_num_TD)), 
                                             max(c(chr_lengths_data$snp_num_TA, 
                                                   chr_lengths_data$snp_num_TD))))
    return(plot)
  } else {
    return(chr_lengths_data)
  }
}

count_SNP_by_window(geno_data = 'source_data/triples_geno_data/2.tsv', plot = T)
count_SNP_by_window(geno_data = 'source_data/triples_geno_data/2.tsv', plot = F)


