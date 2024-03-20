library(dplyr)
library(magrittr)
library(graphics)
data("data_engage")

fit <- try(fn_model(),
           silent = TRUE)

test_that("no error in fitting fn_model for the data", {
  
  expect_false(inherits(fit, "try-error"))
  
})