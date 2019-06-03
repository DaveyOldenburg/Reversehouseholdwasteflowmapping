trashcans_0.1 <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/trashcans_0.1.tsv")
trashcans_0.4 <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/trashcans_0.4.tsv")
trashcans_0.8 <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/trashcans_0.8.tsv")
trashcans_1 <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/trashcans_1.tsv")
View(trashcans_0.4)
View(trashcans_0.1)

trashcans_0.1_bb <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/Bigger buffer/trashcans_2_0.1.tsv")
trashcans_0.4_bb <- read.delim("/Users/daveyoldenburg 1/Desktop/Test for LCA/Validation/Validation Values/total/Bigger buffer/trashcans_2_0.4.tsv")
View(trashcans_0.4_bb)
View(trashcans_0.1_bb)

library(sf)
shape <- read_sf(dsn = "~/Desktop/Thesis/Gemeente Containers/Containerlocations_GM.shp", layer = "Containerlocations_GM")
trashcans_0.1_added <- merge(x = trashcans_0.1, y = shape, by = "LOCATIENR", all.x = TRUE)
trashcans_0.4_added <- merge(x = trashcans_0.4, y = shape, by = "LOCATIENR", all.x = TRUE)
trashcans_0.8_added <- merge(x = trashcans_0.8, y = shape, by = "LOCATIENR", all.x = TRUE)
trashcans_1_added <- merge(x = trashcans_1, y = shape, by = "LOCATIENR", all.x = TRUE)

trashcans_0.1_bb_added <- merge(x = trashcans_0.1_bb, y = shape, by = "LOCATIENR", all.x = TRUE)
trashcans_0.4_bb_added <-merge(x = trashcans_0.4_bb, y = shape, by = "LOCATIENR", all.x = TRUE)

newdata_0.1 <- subset(trashcans_0.1_added, as.numeric(as.character(trashcans_0.1_added$weg_su_res)) > 2000)
newdata_0.4 <- subset(trashcans_0.4_added, as.numeric(as.character(trashcans_0.4_added$weg_su_res)) > 2000)
newdata_0.8 <- subset(trashcans_0.8_added, as.numeric(as.character(trashcans_0.8_added$weg_su_res)) > 2000)
newdata_1 <- subset(trashcans_1_added, as.numeric(as.character(trashcans_1_added$weg_su_res)) > 2000)

View(newdata_0.4)
View(newdata_0.1)

one <- cor(trashcans_0.1_added$count, as.numeric(as.character(trashcans_0.1_added$W1_su_rest)))
two <- cor(trashcans_0.4_added$count, as.numeric(as.character(trashcans_0.4_added$W1_su_rest)))
three <-cor(trashcans_0.8_added$count, as.numeric(as.character(trashcans_0.8_added$W1_su_rest)))
four <- cor(trashcans_1_added$count, as.numeric(as.character(trashcans_1_added$W1_su_rest)))

p <- c(one, two, three, four)
adjustment <- c(0.1,0.4,0.8,1)
install.packages('basicTrendline')
library(basicTrendline)
View(newdata_0.1)
trendline(adjustment, p, model = "line3P", plot = TRUE, linecolor = "blue")

trendline(newdata_0.1$count, as.numeric(as.character(newdata_0.1$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year', linecolor = "red")
trendline(newdata_0.4$count, as.numeric(as.character(newdata_0.4$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year', linecolor = "orange")
trendline(newdata_0.8$count, as.numeric(as.character(newdata_0.8$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year',  linecolor = "green")
trendline(newdata_1$count, as.numeric(as.character(newdata_1$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year', linecolor = "blue")

trendline(trashcans_0.1_bb_added$count, as.numeric(as.character(trashcans_0.1_bb_added$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year', linecolor = "red")
trendline(trashcans_0.4_bb_added$count, as.numeric(as.character(trashcans_0.4_bb_added$weg_su_res)),model = "line2P", plot = TRUE, xlab = 'count of households', ylab = 'kg per year', linecolor = "red")
