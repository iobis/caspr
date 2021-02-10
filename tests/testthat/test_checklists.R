context("checklists")

test_that("the expert priority checklists can be fetched from GitHub", {
  df <- expert_checklists()
  expect_true(is.data.frame(df))
  expect_gt(nrow(df), 0)
})

test_that("the OBIS checklist can be fetched from OBIS", {
  df <- obis_checklist()
  expect_true(is.data.frame(df))
  expect_gt(nrow(df), 0)
})

