mock_vcf = read.csv("inst/annotations.unique.rows.csv", skip = 1)
mock_vcf %>% group_by(Koh476.annotate.class) -> foo
summarize(foo, num89 = n_distinct(Koh89.annotate.class)) |>
  filter(num89 > 1) -> dd
unique(
  dplyr::inner_join(dd, mock_vcf) |>
    select(Koh476.annotate.class, Koh89.annotate.class)
)


split_xx = test_indel_categorization()
split_xx %>% group_by(Koh476.annotate.class) -> foo
summarize(foo, num89 = n_distinct(COSMIC_83)) |> filter(num89 > 1) -> ee

dplyr::inner_join(ee, split_xx) |>
  select(Koh476.annotate.class, COSMIC_83) |>
  distinct() |>
  arrange(Koh476.annotate.class) -> ff


split_xx = test_indel_categorization()
split_xx %>% group_by(Koh89.annotate.class) -> foo
summarize(foo, num89 = n_distinct(COSMIC_83)) |> filter(num89 > 1) -> ee

dplyr::inner_join(ee, split_xx) |>
  select(Koh89.annotate.class, COSMIC_83) |>
  distinct() |>
  arrange(Koh89.annotate.class) -> gg
