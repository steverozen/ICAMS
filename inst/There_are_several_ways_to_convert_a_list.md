There are several ways to convert a list of lists to a tibble. Here are the main alternatives:

  1. Current approach: dplyr::bind_rows() (Already in your code)

  list_of_lists_to_tibble <- function(list_of_lists) {
    dplyr::bind_rows(list_of_lists)
  }
  Pros: Simple, handles missing columns gracefully (fills with NA), returns tibble
  Cons: Requires dplyr dependency

  2. Using purrr::map_df() or purrr::map_dfr()

  list_of_lists_to_tibble <- function(list_of_lists) {
    purrr::map_dfr(list_of_lists, identity)
  }
  # Or more explicitly:
  list_of_lists_to_tibble <- function(list_of_lists) {
    purrr::map_dfr(list_of_lists, ~tibble::as_tibble(.x))
  }
  Pros: Part of tidyverse, functional programming style, map_dfr = "map and bind rows"
  Cons: Requires purrr dependency, essentially does the same as bind_rows() internally

  3. Using base R with do.call(rbind, ...)

  list_of_lists_to_tibble <- function(list_of_lists) {
    tibble::as_tibble(do.call(rbind, lapply(list_of_lists, as.data.frame)))
  }
  Pros: No tidyverse dependencies (except tibble for return type)
  Cons: More verbose, doesn't handle missing columns as gracefully, type coercion issues

  4. Using data.table::rbindlist() (ICAMS already uses data.table)

  list_of_lists_to_tibble <- function(list_of_lists) {
    tibble::as_tibble(data.table::rbindlist(list_of_lists, fill = TRUE))
  }
  Pros: Very fast for large datasets, fill=TRUE handles missing columns, ICAMS already depends on data.table
  Cons: Returns data.table (need to convert to tibble if desired)

  5. Manual construction with tibble::tibble()

  list_of_lists_to_tibble <- function(list_of_lists) {
    # Extract all unique column names
    all_cols <- unique(unlist(lapply(list_of_lists, names)))

    # Create columns by extracting each field
    result <- lapply(all_cols, function(col) {
      sapply(list_of_lists, function(x) x[[col]] %||% NA)
    })
    names(result) <- all_cols

    tibble::as_tibble(result)
  }
  Pros: Full control over the process
  Cons: Much more complex, manual handling of missing values

  Recommendation for ICAMS:

  Since ICAMS already uses data.table extensively (I see it in the code), I'd recommend:

  list_of_lists_to_tibble <- function(list_of_lists) {
    if (length(list_of_lists) == 0) {
      return(data.table::data.table())
    }
    data.table::rbindlist(list_of_lists, fill = TRUE)
  }

  This is:
  - Fast (important for genomics data)
  - No new dependencies (data.table already in ICAMS)
  - Handles missing columns with fill = TRUE
  - Consistent with the rest of ICAMS codebase (returns data.table like line 212)

  The current dplyr::bind_rows() approach works fine, but if you're trying to minimize dependencies or maximize performance, data.table::rbindlist() is a better fit for the ICAMS architecture.      
     