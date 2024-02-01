### A variety of useful functions and helpers for dealing with Synapse/R/python/data weirdness

# get unique vals function

get_unique_vals <- function(data) {
  vals <- purrr::map(data, function(x) {
    unique(x)
  })
  names(vals) <- colnames(data)
  vals
}

# coalesce join function from here: https://alistaire.rbind.io/blog/coalescing-joins/

coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}

# this function takes a synID of the current version of the file you'd like to update, 
# the synID of the file you'd like to upload as the new version of the current file,
# and an optional parameter to add to the version comment

upload_new_version <- function(current_id, new_id, version_comment = NULL) {
  # download current and new files
  current_obj <- synGet(current_id)
  new_obj <- synGet(new_id)
  
  # edit current file path to point to new version
  current_obj$path <- new_obj$path
  
  
  current_obj$properties$versionComment <- version_comment
  synStore(current_obj)
}

# remove JSON from multi-annotations
#' @description Remove the JSON-specific elements from a string. For example,
#' a JSON string could look like `"[\"foo\", \"bar\"]"`. This function would
#' return the string `"foo, bar"`.
#'
#' @param json_string Character string with JSON elements.
#' @param remove_spaces `TRUE` to also remove spaces, else `FALSE`. Note that
#' this will remove *all* spaces. Default is `FALSE`.
#' @return Character string without brackets, escaped quotes, and (optionally)
#' spaces.

clean_json_string <- function(json_string, remove_spaces = FALSE) {
  if (remove_spaces) {
    return(gsub("\\[|\"|\\]| ", "", json_string))
  } else {
    return(gsub("\\[|\"|\\]|", "", json_string))
  }
}

# format a JSON string for multi-value annotations
#' @description ADD JSON-specific elements to a string to format for multi-value annotations set via a file view. For example,
#' providing the string `"foo, bar"` should return the JSON string `"[\"foo\", \"bar\"]"`. 
#'
#' @param string Character string.
#' @param has_spaces `TRUE` if the provided string has spaces between the elements, otherwise `FALSE`. Default is `FALSE`.
#' @return Character string with brackets, escaped quotes, and (optionally)
#' spaces.
make_json_string <- function(string, has_spaces = TRUE) {
  ends <- paste0("[\"", string, "\"]")
  if (has_spaces){
    new_json_string <- gsub(", ", "\\\", \\\"", ends)
  } else {
    new_json_string <- gsub(",", "\\\",\\\"", ends)
  }
  return(new_json_string)
}

### This version of the make_json_string function uses purrr to handle 
### NA values -- NAs are untouched, strings are converted to json
make_json_string_purrr <-
  function(string) {
    purrr::map_chr(string, function(x) {
      if (is.na(x)) {
        NA_character_
      } else {
        make_json_string(x)
      }
    })
  }

# Many ROSMAP projids have 1 or 2 leading zeros, which tend to get dropped 
# when downloading annotations or metadata csvs. This function takes as input
# a dataframe or tibble with a column labeled `projid` that contains ROSMAP projids
# as character values, appends the appropriate number of leading 0s to any projids 
# with fewer than 8 characters, and returns the original dataframe with the fixed
# projid values. 

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

