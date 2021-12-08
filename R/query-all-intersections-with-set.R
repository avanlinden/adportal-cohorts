# function to apply upset queries that select all intersections with one ore more terms

# test plot

upset(
  rosmap_boolean_categories,
  rev(sortedUpsetCategories),
  min_size = 30,
)

# all possible intersections given data, categories, and minimum intersection size
intersections <-  unique(upset_data(data = rosmap_boolean_categories,
                                    intersect = sortedUpsetCategories,
                                    min_size = 30)$plot_intersections_subset)


# function that takes the data, sets, and minimum intersection size for the upset plot data,
# plus a vector of strings matching one or more sets, and returns a vector of every intersection
# for the given upset plot that contains at least one of those sets

define_query_intersections <- function(data, sets, min_size = 0,  sets_to_query) {
  intersections <- unique(upset_data(data, sets, min_size)$plot_intersections_subset)
  intersection_query_by_set <- paste(sets_to_query, collapse = "|")
  queried_intersections <- intersections[str_detect(intersections, intersection_query_by_set)]
  return(queried_intersections)
}

test <- define_query_intersections(rosmap_boolean_categories, 
                                   sortedUpsetCategories, 
                                   min_size = 30,
                                   sets_to_query = c("peripheral metabolomics", "monocyte bulk RNAseq"))


# make a list out of each of these queries


queries <-  list(
  upset_query(
    intersect = c("brain bulk RNAseq", "genomic variants"),
    color = 'red',
    fill = 'red'
  ),
  upset_query(
    set = 'genomic variants',
    fill = 'blue'
  )
)
queries

test_list<- map(test, ~list(set = NULL,
                            intersect = unlist(str_split(.x, "-")), 
                            group_by_group = NULL,
                            only_components = NULL,
                            color = "red",
                            fill = "red"))


glue_list <- map_depth(test_list, 0, ~glue::glue("upset_query({.x})"))

obj <- upset_data(rosmap_boolean_categories, sortedUpsetCategories, min_size = 20)

test_intersections <- map(test, ~list(glue::glue("upset_query({intersect = unlist(str_split(.x, \"-\"))})")))

test_glue <- glue::glue("upset_query({test_intersections}")

list_glue_list <- as.list(glue_list)

upset_query_list <- map(test, ~upset_query(intersect = str_replace(unlist(str_split(.x, "-")), "_", "-"),
                                           fill = "purple"))

str(obj$matrix_frame$intersection)

str(obj$matrix)

try_one <- upset_query_list[1]
try_one

try_queries <- map(upset_query_list, .x)

fix <- str_replace(test_intersect[[1]], "_", "-")

set_queries <- list(upset_query(set = "LC-SRM", fill = "purple"),
                    upset_query(set = "TMT quantitation", fill = "purple"))

all_queries <- append(upset_query_list, set_queries)
