##### AUTHOR:     JILL THORR, M.S.
##### SCRIPT: 		function_loadRData.R
##### DATE:		  10 DEC 2022
##### OBJECTIVE:	Load R data and assign any name to the 
#####				associated file

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
