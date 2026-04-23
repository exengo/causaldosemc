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

test_that("washout treats unobserved (NA) dose cells as potentially active", {
  dose_matrix <- rbind(
    c(0, NA, 0, 0, 0),
    c(0, 0, 0, NA, 0)
  )

  mask <- causaldosemc:::cdmc_build_eligible_mask(
    dose_matrix = dose_matrix,
    zero_tolerance = 1e-8,
    washout = 2
  )

  # Row 1: cell at t=2 is NA; subsequent observed zero-dose cells within
  # washout window must be excluded because we can't rule out treatment.
  expect_true(mask[1, 1])
  expect_false(mask[1, 2])  # NA itself is never eligible
  expect_false(mask[1, 3])
  expect_false(mask[1, 4])
  expect_true(mask[1, 5])

  # Row 2: NA at t=4; t=5 within washout window must be excluded.
  expect_true(mask[2, 1])
  expect_true(mask[2, 2])
  expect_true(mask[2, 3])
  expect_false(mask[2, 4])  # NA itself is never eligible
  expect_false(mask[2, 5])
})
