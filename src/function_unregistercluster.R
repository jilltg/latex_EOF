## SCRIPT: 		function_unregistercluster.R
## AUTHOR:		JILL THORR, M.S. (THOMPSON-GRIM)
## DATE:   		4 APRIL 2023
## OBJECTIVE:	UNREGISTER THE CLUSTER CREATED WITH
##			     	doParallel::registerDoParallel(cores = detectCores())
##				    THIS CAN BE RUN AFTER FOREACH LOOPS TO 
##				    PREVENT ZOMBIE TASKS

unregister <- function() {
  library(foreach)
  library(doParallel)
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

