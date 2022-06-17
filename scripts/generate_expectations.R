library(straf)

input_diploid <- list(datapath="./tests/input_files/exampleSTRAFdiplo.txt")
dat_diploid <- straf::createGenind(input_diploid, ploidy = 2)
freq_diploid <- straf::getFreqAllPop(dat_diploid)
indices_diploid <- straf::getIndicesAllPop(dat_diploid)


input_haploid <- list(datapath="./tests/input_files/exampleSTRAFhaplo.txt")
dat_haploid <- straf::createGenind(input_haploid, ploidy = 1)
freq_haploid <- straf::getFreqAllPop(dat_diploid)
indices_haploid <- straf::getIndicesAllPop(dat_diploid)


input_haploid_point <- list(datapath="./tests/input_files/exampleSTRAFhaplo_point.txt")
dat_haploid_point <- straf::createGenind(input_haploid_point, ploidy = 1)
freq_haploid_point <- straf::getFreqAllPop(dat_diploid)
indices_haploid_point <- straf::getIndicesAllPop(dat_diploid)

expectations <- list(
  dat_diploid=dat_diploid,
  freq_diploid=freq_diploid,
  indices_diploid=indices_diploid,
  dat_haploid=dat_haploid,
  freq_haploid=freq_haploid,
  indices_haploid=indices_haploid,
  dat_haploid_point=dat_haploid_point,
  freq_haploid_point=freq_haploid_point,
  indices_haploid_point=indices_haploid_point
)

saveRDS(expectations, file = "./tests/expectations.rds")