require(stringi)
require(stringr)
require(data.table)
require(dplyr)

setwd("~/Documents/cenrich/")

## glm (N ~ Germline Variants + PC1 + PC2 + PC3 + PC4, family= ''binomial'')
## where: N = case (1) or control (0), Germline Variants = indicating the number of samples that carries rare pathogenic germline variants for each gene-tissue pair, PC values from PCA analysis of PCAWG and 1KG to control population structures were used as an input for the regression. 

inputD <- fread(inputdata, data.table=F)
unique(inputD$hist)
nrow(inputD) 

calist <- "Pancan"

RGLM <- function(caa){
  
  print(caa)
  
  inputD_c <- inputD[which(inputD$hist == caa),]
  inputD_n <- inputD[which(inputD$hist == "TG"),]
  inputR <- inputD
  
  case_add <- as.data.frame(matrix(c("case-fake-sample", 1, caa, mean(inputR$PC1), mean(inputR$PC2), mean(inputR$PC3), mean(inputR$PC4), 
                                     rep(1, ncol(inputR)-7)), nrow=1), stringsAsFactors=F)
  control_add <- as.data.frame(matrix(c("control-fake-sample", 0, "TG", mean(inputR$PC1), mean(inputR$PC2), mean(inputR$PC3), mean(inputR$PC4),
                                        rep(1, ncol(inputR)-7)), nrow=1), stringsAsFactors=F)
  names(case_add) <- colnames(inputR)
  names(control_add) <- colnames(inputR)
  F.input <- rbind(inputR, case_add, control_add)
  print(nrow(F.input))
  
  F.input$type1 <- as.numeric(F.input$type1)
  F.input[8:ncol(F.input)] <- lapply(F.input[8:ncol(F.input)], as.integer)
  F.input$PC1 <- as.numeric(F.input$PC1)
  F.input$PC2 <- as.numeric(F.input$PC2)
  F.input$PC3 <- as.numeric(F.input$PC3)
  F.input$PC4 <- as.numeric(F.input$PC4)
  
  vglist <- names(F.input)[8:ncol(F.input)]
  
  glmmm <- lapply(vglist, function(x) {
    
    { glm.codee <- glm(substitute(type1 ~ i + PC1 + PC2 + PC3 + PC4, list(i = as.name(x))), data=F.input, family='binomial')
    
    ###variable ( 'i'= variant carrier / LOH / Disease / Pathway )
    
    }})
  
  glmmm.s <- lapply(glmmm, summary)
  glmmm.s_r <- lapply(glmmm.s, '[[', 'coefficients')
  glmmm.s_r <- mapply(cbind, glmmm.s_r, SIMPLIFY=F)
  glmmm.s_tot <- do.call(rbind, glmmm.s_r)
  
  return(glmmm.s_tot)    
}

caEnrich <- lapply(calist, RGLM) %>% as.data.frame()
caEnrich
caEnrich$OR <- exp(as.numeric(caEnrich$Estimate))
caEnrich$FDR <- p.adjust(caEnrich$Pr...z.., method = 'fdr')
head(caEnrich)

write.table(caEnrich, "Cancer_Enrich.sig.can.table.tsv", sep = "\t")
