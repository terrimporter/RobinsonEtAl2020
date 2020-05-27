#Chloe Robinson, Feb 28 2020

install.packages('VennDiagram')
install.packages('RAM')
library(VennDiagram)
install.packages("venneuler")
library(venneuler)

#Use values from datafile 'Values_for_VennDiagrams.xlsx'


#Single sampling + Treatment (Ethanol vs Antifreeze)

grid.newpage()
draw.pairwise.venn(area1 = 2+22, area2 = 4+22, cross.area = 22, category = c("Ethanol", "Antifreeze"), lty = "blank", label = TRUE, cex = 2, cat.cex = 1.5, fill = c("light blue", "tomato"))

#Single sampling + Site (Site 1 vs Site 2 vs Site 3)

grid.newpage()
draw.triple.venn(area1 = 3+1+14+2, area2 = 1+3+14+1, area3= 4+2+14+3, n12 = 14+1, n23 = 3+14, n13 = 2+14, n123 = 14, category = c("Beaver 18", "Beaver 19", "Clair 12"), lty = "blank",cex = 2, cat.cex = 1.5, fill = c("gold", "darkorchid1", "light green"))

#Single sampling + Method (Antifreeze only; Method 1(no evap) vs Method 2 (evap))

grid.newpage()
draw.pairwise.venn(area1 = 4+23, area2 = 1+23, cross.area = 23, category = c("Method 1","Method 2"), lty = "blank", cex = 2, cat.cex = 1.5,fill = c("lightpink", "chocolate1"))

#Paired sampling + Treatment (Ethanol vs Antifreeze)

grid.newpage()
draw.pairwise.venn(area1 = 4+26, area2 = 5+26, cross.area = 26, category = c("Ethanol", "Antifreeze"), lty = "blank", label = TRUE, cex = 2, cat.cex = 1.5, fill = c("light blue", "tomato"))

#Paired sampling + Site (Site 1 vs Site 2 vs Site 3)

grid.newpage()
draw.triple.venn(area1 = 2+5+13, area2 = 5+7+13+5, area3= 3+7+13, n12 = 5+13, n23 = 7+13, n13 = 13, n123 = 13, category = c("Laurel 4", "Laurel 7", "Laurel 10"), lty = "blank",cex = 2, cat.cex = 1.5, fill = c("gold", "darkorchid1", "light green"))

#Paired sampling + Method (Antifreeze only; Method 1(no evap) vs Method 2 (evap))

grid.newpage()
draw.pairwise.venn(area1 = 1+27, area2 = 4+27, cross.area = 27, category = c("Method 1","Method 2"), lty = "blank", cex = 2, cat.cex = 1.5,fill = c("lightpink", "chocolate1"))

