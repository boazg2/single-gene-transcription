library(ggplot2)
library(dplyr)

data = data.frame(
  gene = c(
    "lacPs", "lacPs", "lacPs", "lacPs", "lacPs", "lacPs", "lacPs",
    "lacPsd", "lacPsd", "lacPsd", "lacPsd", "lacPsd", "lacPsd", "lacPsd", "lacPsd", "lacPsd",
    "tyrT", "tyrT", "tyrT", "tyrT", "tyrT", "tyrT",
    "tyrTd", "tyrTd", "tyrTd", "tyrTd", "tyrTd", "tyrTd",
    "galP", "galP", "galP", "galP", "galP", "galP", "galP", "galP", "galP"
  ),
  sigma = c(0,.016,.032,.042,.057,.062,.081,
            0,.012,.023,.035,.045,.056,.069,.081,.092,
            0,.015,.024,.036,.054,.068,
            0,.015,.024,.036,.054,.068,
            0,.015,.024,.038,.05,.059,.071,.079,.099),
  expr = c(1,2,7,12,19,22,13,
               8,10,20,33,36,32,26,23,21,
               0,0,4,10,140,100,
               0,15,50,120,95,100,
               1,2.3,4.1,6.4,11,12.5,11.5,9.7,8)
) %>%
  group_by(gene) %>%
  mutate(sigma = -sigma, rel_expr = expr / max(expr)) %>%
  ungroup

ggplot(data = data, aes(x = sigma, y = rel_expr, color = gene)) + 
  geom_line() +
  geom_point() +
  theme_classic() +
  ylab("Relative Expression") + 
  xlab("Sigma")
