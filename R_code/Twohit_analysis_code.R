require(stringi)
require(stringr)
require(data.table)
require(dplyr)

setwd("~/Documents/cenrich/")

## glm (N ~ Germline Variants + PC1 + PC2 + PC3 + PC4, family= ''binomial'')
## where: N represents samples with LOH event (1) or samples without LOH event (0), Germline Variants indicates the number of samples carrying rare pathogenic variants (Tier1 + Tier2) for each gene-tissue pair, PC values from PCA analysis of PCAWG only were used as an input for the regression. 

inputL <- fread(inputloh, data.table=F)
inputL$histology_abbreviation <- "Pancan"

glmmm_T <- function(catype){
  
  print(catype)
  
  inputL2 <- inputL[which(inputL$histology_abbreviation == catype),]
  nrow(inputL2)
  case_add <- as.data.frame(matrix(c("case-fake-sample", 1, catype, mean(inputL2$PCAWG_PC1), mean(inputL2$PCAWG_PC2), mean(inputL2$PCAWG_PC3), mean(inputL2$PCAWG_PC4),
                                     rep(1, ncol(inputL)-7)), nrow=1), stringsAsFactors=F)
  control_add <- as.data.frame(matrix(c("control-fake-sample", 1, catype, mean(inputL2$PCAWG_PC1), mean(inputL2$PCAWG_PC2), mean(inputL2$PCAWG_PC3), mean(inputL2$PCAWG_PC4), 
                                        rep(0, 4), rep(1, 4)), nrow=1), stringsAsFactors=F)
  
  ncol(inputL2)
  
  names(case_add) <- colnames(inputL2)
  names(control_add) <- colnames(inputL2)
  inputL2 <- rbind(inputL2, case_add, control_add)
  print(nrow(inputL2))
  
  inputL2$`Case-Control` <- as.numeric(inputL2$`Case-Control`)
  inputL2[8:ncol(inputL2)] <- lapply(inputL2[8:ncol(inputL2)], as.integer)
  inputL2$PCAWG_PC1 <- as.numeric(inputL2$PCAWG_PC1)
  inputL2$PCAWG_PC2 <- as.numeric(inputL2$PCAWG_PC2)
  inputL2$PCAWG_PC3 <- as.numeric(inputL2$PCAWG_PC3)
  inputL2$PCAWG_PC4 <- as.numeric(inputL2$PCAWG_PC4)

  varlist <- colnames(inputL2)[12:15]
  print(length(varlist))
  
  glmmm_T <- lapply(varlist, function(x) {
    suppressWarnings({
      
      s = x
      
      formul <- paste0(s,"_LOH"," ~ ",x, " + PCAWG_PC1 + PCAWG_PC2 + PCAWG_PC3 + PCAWG_PC4")
      glmmm <- glm(as.formula(formul), data=inputL2, family='binomial')
      
    })})
  
  glmmm.s <- lapply(glmmm_T, summary)
  glmmm.s_r <- lapply(glmmm.s, '[[', 'coefficients')
  glmmm.s_r <- mapply(cbind, glmmm.s_r, SIMPLIFY=F)
  caEnrichT <- do.call(rbind, glmmm.s_r) %>% as.data.frame()
  
  return(caEnrichT)    
}

TEnrich <- lapply(inputL$histology_abbreviation, glmmm_T) %>% as.data.frame()
TEnrich
TEnrich$OR <- exp(as.numeric(TEnrich$Estimate))
TEnrich$FDR <- p.adjust(TEnrich$Pr...z.., method = 'fdr')
