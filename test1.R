library(ggplot2)

#library(pestim)
profvis({

x <- seq(0.10, 0.65, by = 0.01)
y <- dtriang(x, xMin = 0.10, xMax = 0.65, xMode = 0.32)
df <- data.frame(x = x, y= y)
ggplot(df, aes(x, y)) + geom_line() + scale_x_continuous(limits = c(0, 1)) + xlab('u') + ylab('Probability Density')
})
