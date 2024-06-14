balance_ratio <- function(x){
  table(x) -> counts
  sort(counts) |> rev() -> sorted_counts
  sorted_counts[1]/sorted_counts[length(sorted_counts)]
}
