test_that("washout masking drops immediate post-treatment zero-dose cells", {
  dose_matrix <- rbind(
    c(0, 1, 0, 0, 0),
    c(0, 0, 0, 0, 0)
  )

  mask <- causaldosemc:::cdmc_build_eligible_mask(
    dose_matrix = dose_matrix,
    zero_tolerance = 1e-8,
    washout = 2
  )

  expect_true(mask[1, 1])
  expect_false(mask[1, 3])
  expect_false(mask[1, 4])
  expect_true(mask[1, 5])
})
