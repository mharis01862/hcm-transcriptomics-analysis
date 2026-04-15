library(ggplot2)
library(readr)
library(dplyr)
library(forcats)

david_data <- read_csv("david_chart.csv")
colnames(david_data) <- make.names(colnames(david_data))

david_data <- david_data %>%
  mutate(NegLog10PValue = -log10(PValue)) %>%
  arrange(NegLog10PValue) %>%
  mutate(Term = factor(Term, levels = Term))

category_colors <- c(
  "GOTERM_BP_DIRECT" = "#3ef0bd",
  "GOTERM_CC_DIRECT" = "#ebb48a",
  "GOTERM_MF_DIRECT" = "#d1f598"
)

ggplot(david_data, aes(x = NegLog10PValue,
                       y = fct_rev(Term),
                       fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(title = "GO Functional Analysis",
       x = "-log10(PValue)")
