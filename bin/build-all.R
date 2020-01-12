#!/usr/bin/env Rscript

here <- getwd()
root <- dirname(here)
mydir <- file.path(root,"analysis")
scripts <- list.files(mydir,pattern="build",full.names=TRUE)
for (script in scripts){
	system(script,ignore.stdout=FALSE,ignore.stderr=FALSE)
}
