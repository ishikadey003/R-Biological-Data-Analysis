# Load the PlantGrowth dataset
data("PlantGrowth")
# View first few rows
head(PlantGrowth)
# Check structure
str(PlantGrowth)
View(PlantGrowth)
# Mean weight of each group
mean_ctrl <- mean(PlantGrowth$weight[PlantGrowth$group=="ctrl"])
mean_trt1 <- mean(PlantGrowth$weight[PlantGrowth$group=="trt1"])
mean_trt2 <- mean(PlantGrowth$weight[PlantGrowth$group=="trt2"])
mean_ctrl
mean_trt1
mean_trt2
boxplot(weight ~ group,
data=PlantGrowth,
main="Plant Growth",
xlab="Treatment",
ylab="Weight")
hist(PlantGrowth$weight,
main="Histogram of Plant Weight",
xlab="Weight")
hist(PlantGrowth$weight,
main = "Plant Growth",
xlab = "Weight",
ylab = "Number of Plants")
library(vioplot)
vioplot(
PlantGrowth$weight[PlantGrowth$group=="ctrl"],
PlantGrowth$weight[PlantGrowth$group=="trt1"],
PlantGrowth$weight[PlantGrowth$group=="trt2"],
names=c("Control","Treatment1","Treatment2")
)
result <- aov(weight ~ group,
data=PlantGrowth)
summary(result)

