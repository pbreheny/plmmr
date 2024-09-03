create_design <- function(dat_file,
                          rds_dir,
                          ...) {
  obj <- readRDS(dat_file)
  switch (class(obj),
    processed_plink = create_design_plink(dat_file = dat_file, rds_dir = rds_dir, ...),
    processed_delim = create_design_delim(dat_file = dat_file, rds_dir = rds_dir, ...),
    processed_matrix = create_design_matrix(dat_file = dat_file, rds_dir = rds_dir, ...),
  )
}