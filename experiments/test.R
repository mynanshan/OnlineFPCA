x <- rnorm(100)
mean(x)

library(ggplot2)

df <- data.frame(
  x = rnorm(100),
  y = rnorm(100)
)

ggplot(df, aes(x, y)) +
  geom_point()
