require(stringi)
require(stringr)
require(data.table)
require(dplyr)

setwd("~/Documents/cenrich/")

## glm (N ~ Disease classes + PC1 + PC2 + PC3 + PC4, family= ''binomial'')
## where: N = case (1) or control (0), Disease classes = indicating the number of samples carrying a rare pathogenic variant for each disease-related gene, PC values from PCA analysis of PCAWG and 1KG to control population structures were used as an input for the regression. 

inputDisease <- fread(inputdatadisease, data.table=F)
unique(inputDisease$hist)
nrow(inputDisease) 

inputDisease$histology_abbreviation <- "Pancan"

controlinput <- inputDisease[which(inputDisease$histology_abbreviation == "1KG"),]
caseinput <- subset(inputDisease, histology_abbreviation != "1KG")
controlinput <- subset(inputDisease, histology_abbreviation == "1KG")

case_add <- as.data.frame(matrix(c("case-fake-sample", 1, "Pancan", mean(inputDisease$PC1), mean(inputDisease$PC2), 
                                   mean(inputDisease$PC3), mean(inputDisease$PC4), rep(1, ncol(inputDisease)-7)), nrow=1), stringsAsFactors=F)
control_add <- as.data.frame(matrix(c("control-fake-sample", 0, "1KG", mean(inputDisease$PC1), mean(inputDisease$PC2), 
                                      mean(inputDisease$PC3), mean(inputDisease$PC4), rep(1, ncol(inputDisease)-7)), nrow=1), stringsAsFactors=F)
names(case_add) <- colnames(inputDisease)
names(control_add) <- colnames(inputDisease)
DF.input <- rbind(inputDisease, case_add, control_add)

DF.input$`Case-Control` <- as.numeric(DF.input$`Case-Control`)
DF.input[8:ncol(DF.input)] <- lapply(DF.input[8:ncol(DF.input)], as.integer)
DF.input$PC1 <- as.numeric(DF.input$PC1)
DF.input$PC2 <- as.numeric(DF.input$PC2)
DF.input$PC3 <- as.numeric(DF.input$PC3)
DF.input$PC4 <- as.numeric(DF.input$PC4)

varlist <- names(DF.input)[8:ncol(DF.input)]

glmmm_D <- lapply(varlist, function(x) {
  suppressWarnings({
    
    glmm <- glm(substitute(`Case-Control` ~ i + PC1 + PC2 + PC3 + PC4   , list(i = as.name(x))), data=DF.input, family='binomial')
    
  })})

glmmm.s <- lapply(glmmm_D, summary)
glmmm.s_r <- lapply(glmmm.s, '[[', 'coefficients')
glmmm.s_r <- mapply(cbind, glmmm.s_r, SIMPLIFY=F)
caEnrichD <- do.call(rbind, glmmm.s_r) %>% as.data.frame()

caEnrichD$OR <- exp(as.numeric(caEnrichD$Estimate))
caEnrichD$FDR <- p.adjust(caEnrichD$Pr...z.., method = 'fdr')
head(caEnrichD)

write.table(caEnrichD, "Disease_analysis.sig.table.tsv", sep = "\t")
