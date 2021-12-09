# Fix leading zeros in ROSMAP projids

fix_leading_zeros <- function(df) {
  fixed_projids <- purrr::map(df$projid, function(x) {
    if (nchar(x) < 8) {
      zeros <- c(glue::glue_collapse(rep("0", (8 - nchar(x)))))
      x <- paste0(zeros, x)
    }
    return(x)
  }
  )
  unnest(mutate(df, projid = fixed_projids), projid)
}
