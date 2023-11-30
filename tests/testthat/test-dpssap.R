test_that("Test for dpssap", {
  dw <- multitaper::dpss(n = 100, k = 12, nw = 6)$v
  dwap <- dpssap(V = dw, maxdeg = 3)
  expect_equal(class(dwap), "list")
  expect_length(dwap, 3)
  expect_error(dpssap(V=dw))
})
