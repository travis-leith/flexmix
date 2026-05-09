test_that("3-component mixture with 2 concomitant variables recovers true coefficients", {
  set.seed(42)
  n <- 2000

  # Concomitant variables
  z1 <- rnorm(n)
  z2 <- rnorm(n)

  # True concomitant coefficients (for multinomial logit, component 1 is reference)
  # Component 2 vs 1: intercept=-0.5, z1=1.5, z2=-1.0
  # Component 3 vs 1: intercept=0.5,  z1=-1.0, z2=1.5
  log_p2 <- -0.5 + 1.5 * z1 - 1.0 * z2
  log_p3 <- 0.5 - 1.0 * z1 + 1.5 * z2

  denom <- 1 + exp(log_p2) + exp(log_p3)
  p1 <- 1 / denom
  p2 <- exp(log_p2) / denom
  p3 <- exp(log_p3) / denom

  # Assign components
  component <- apply(cbind(p1, p2, p3), 1, function(p) sample(1:3, 1, prob = p))

  # True regression coefficients per component
  # y = intercept + beta*x + noise
  true_intercepts <- c(0.0, 5.0, 10.0)
  true_slopes <- c(2.0, -1.5, 3.0)
  true_sigma <- 1.0

  x <- rnorm(n)
  y <- true_intercepts[component] +
    true_slopes[component] * x +
    rnorm(n, sd = true_sigma)

  dat <- data.frame(y = y, x = x, z1 = z1, z2 = z2)

  # Fit model
  m <- flexmix_fast(
    y ~ x,
    data = dat,
    k = 3,
    concomitant = FLXPmultinom(~ z1 + z2),
  )


  # Extract component regression coefficients (intercept + slope per component)
  comp_coefs <- sapply(1:3, function(k) {
    m@components[[k]][[1]]@parameters$coef
  })
  # comp_coefs is 2 x 3: rows = (intercept, slope), cols = components

  # Match estimated components to true components by intercept ordering
  est_intercepts <- comp_coefs[1, ]
  order_est <- order(est_intercepts)
  order_true <- order(true_intercepts)

  matched_intercepts <- est_intercepts[order_est]
  matched_slopes <- comp_coefs[2, order_est]

  expect_equal(
    matched_intercepts,
    true_intercepts[order_true],
    tolerance = 0.3,
    label = "Estimated intercepts should be close to true intercepts"
  )
  expect_equal(
    matched_slopes,
    true_slopes[order_true],
    tolerance = 0.3,
    label = "Estimated slopes should be close to true slopes"
  )

  
})
