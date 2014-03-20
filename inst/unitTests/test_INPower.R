test_INPower <- function() {

  MAFs  <- seq(0.1, 0.5, by=0.01)
  betas <- seq(0.05, 0.45, by=0.01) 
  pow   <- seq(0.5, 0.9, by=0.01)
  vec1  <- c(24.9, 34.1, 38.7, 41.6, 43.8)
  vec2  <- c(2.9, 13.3, 20.4, 24.9, 28.1)

  res   <- INPower(MAFs, betas, pow, span=0.5, binary.outcome=FALSE, 
           sample.size=seq(1000,5000,by=1000), 
           signif.lvl=10^(-7), k=seq(5, 20, by=5))
  checkEqualsNumeric(res$future.study.summary$e.discov[, 2], vec1, tolerance=1.0e-4)

  res   <- INPower(MAFs, betas, pow, span=0.5, binary.outcome=TRUE, 
           sample.size=seq(1000,5000,by=1000), 
           signif.lvl=10^(-7), k=seq(5, 20, by=5))
  checkEqualsNumeric(res$future.study.summary$e.discov[, 2], vec2, tolerance=1.0e-4)


}
