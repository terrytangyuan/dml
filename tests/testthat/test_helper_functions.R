context('helper functions')

test_that('package set up successfully', {
  expect_that(sum(1,2), not(throws_error()))
})
