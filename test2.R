library(pestim)
profvis({
alphas <- c(1, 5, 10, 100, 1000)

mode <- 35
df <- lapply(alphas, function(alpha){

  x <- 0:100
  y <- dgamma(x, shape = alpha + 1, scale = mode / alpha)
  z <- as.character(alpha)
  output <- data.frame(x = x, y = y, alpha = z)
  return(output)
})
df <- Reduce(rbind, df)
ggplot(df, aes(x, y, col = alpha, group = alpha)) + geom_line(aes(linetype = alpha)) + scale_x_continuous(limits = c(0, 100)) + xlab('') + ylab('')
})
