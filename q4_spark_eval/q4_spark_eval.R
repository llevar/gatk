install.packages("ggplot2")
library(ggplot2)

test_case_2 <- read.csv("/Users/tom/workspace/gatk/test_case_2_results.csv")
test_case_2$Total.cores <- with(test_case_2, Number.of.executors * Executor.cores)
test_case_2$cm <- with(test_case_2, Total.cores / Time..mins.)
test_case_2$Executor.cores <- as.factor(test_case_2$Executor.cores)

png("/Users/tom/workspace/gatk/test_case_2_time.png")
qplot(Total.cores, Time..mins., data = test_case_2, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Time (mins)",
      main = "Time to run count reads (test case 2) by total number of cores")
dev.off()

install.packages("dplyr")
library("dplyr")

speedup <- test_case_2 %>%
  filter(Total.cores == 4) %>%
  summarize(mean4cores = mean(Time..mins.)) %>%
  as.numeric()
test_case_2$Speedup <- with(test_case_2, 4 * speedup / Time..mins.)
png("/Users/tom/workspace/gatk/test_case_2_speedup.png")
qplot(Total.cores, Speedup, data = test_case_2, geom = c("point", "line"), colour = Executor.cores,
      xlab = "Cores", ylab = "Speedup",
      main = "Count reads (test case 2) speedup by total number of cores")
dev.off()

