library(straf)
library(testthat)

input_diploid <- list(datapath="./exampleSTRAFdiplo.txt")
dat_diploid <- straf::createGenind(input_diploid, ploidy = 2)
freq_diploid <- straf::getFreqAllPop(dat_diploid)
indices_diploid <- straf::getIndicesAllPop(dat_diploid)

input_diploid2 <- list(datapath="./exampleSTRAFdiplo_empty_lines.txt")
dat_diploid2 <- straf::createGenind(input_diploid2, ploidy = 2)
freq_diploid2 <- straf::getFreqAllPop(dat_diploid2)
indices_diploid2 <- straf::getIndicesAllPop(dat_diploid2)

input_haploid <- list(datapath="./exampleSTRAFhaplo.txt")
dat_haploid <- straf::createGenind(input_haploid, ploidy = 1)
freq_haploid <- straf::getFreqAllPop(dat_haploid)
indices_haploid <- straf::getIndicesAllPop(dat_haploid)

input_haploid_point <- list(datapath="./exampleSTRAFhaplo_point.txt")
dat_haploid_point <- straf::createGenind(input_haploid_point, ploidy = 1)
freq_haploid_point <- straf::getFreqAllPop(dat_haploid_point)
indices_haploid_point <- straf::getIndicesAllPop(dat_haploid_point)

expected <- readRDS("./expectations.rds")

test_that("STRAF file parser", {
  expect_true(compare(dat_diploid, expected$dat_diploid)$equal)
  expect_true(compare(dat_diploid2, expected$dat_diploid2)$equal)
  expect_true(compare(dat_haploid, expected$dat_haploid)$equal)
  expect_true(compare(dat_haploid_point, expected$dat_haploid_point)$equal)
})

test_that("STRAF allele frequencies computation", {
  expect_true(compare(freq_diploid, expected$freq_diploid)$equal)
  expect_true(compare(freq_diploid2, expected$freq_diploid2)$equal)
  expect_true(compare(freq_haploid, expected$freq_haploid)$equal)
  expect_true(compare(freq_haploid_point, expected$freq_haploid_point)$equal)
})

test_that("STRAF indices computation", {
  expect_true(compare(indices_diploid, expected$indices_diploid)$equal)
  expect_true(compare(indices_diploid2, expected$indices_diploid2)$equal)
  expect_true(compare(indices_haploid, expected$indices_haploid)$equal)
  expect_true(compare(indices_haploid_point, expected$indices_haploid_point)$equal)
})
