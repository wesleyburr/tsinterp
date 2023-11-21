test_that("Generate ttp and freqsave test", {
  data("flux")
  z1 <- flux$SagOrig
  z1 <- z1[1:50] 
  gen_result <- gen_m_t(zI = z1, N=length(z1), delT=10, sigClip = 0.05, progress= FALSE)
  expect_length(gen_result, 2)
  expect_named(gen_result, c("TtP", "freqSave"))
  expect_equal(names(gen_result[1]), "TtP")
})
