.libPaths(c('C:/Code/UG_ERP_Project/renv/library/R-4.3/x86_64-w64-mingw32',.libPaths()))
library(glmmTMB);library(emmeans)
f<-getFromNamespace('test.emmGrid','emmeans');cat(paste(deparse(f),collapse='\n'))
cat('\nHELPERS\n');print(grep('test',ls(getNamespace('emmeans'),all.names=TRUE),value=TRUE))
