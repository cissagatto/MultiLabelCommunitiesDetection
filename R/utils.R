##################################################################################################
# Tests Communities Detection Methods                                                            #
# Copyright (C) 2021                                                                             #
#                                                                                                #
# This code is free software: you can redistribute it and/or modify it under the terms of the    #
# GNU General Public License as published by the Free Software Foundation, either version 3 of   #
# the License, or (at your option) any later version. This code is distributed in the hope       #
# that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    #
# more details.                                                                                  #
#                                                                                                #
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin                     #
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus Sao Carlos           #
# Computer Department (DC: https://site.dc.ufscar.br/)                                           #
# Program of Post Graduation in Computer Science (PPG-CC: http://ppgcc.dc.ufscar.br/)            #
# Bioinformatics and Machine Learning Group (BIOMAL: http://www.biomal.ufscar.br/)               #
#                                                                                                #
##################################################################################################


##################################################################################################
# Configures the workspace according to the operating system                                     #
##################################################################################################
sistema = c(Sys.info())
FolderRoot = ""
if (sistema[1] == "Linux"){
  FolderRoot = paste("/home/", sistema[7], "/TestsCommunitiesDetection", sep="")
  setwd(FolderRoot)
} else {
  FolderRoot = paste("C:/Users/", sistema[7], "/TestsCommunitiesDetection", sep="")
  setwd(FolderRoot)
}
setwd(FolderRoot)

##################################################################################################

folders <- function(dataset_name){
  retorno = list()
  
  FolderScripts =  paste(FolderRoot, "/R", sep="")
  if(dir.exists(FolderScripts)==FALSE){dir.create(FolderScripts)}
  
  FolderDocs =  paste(FolderRoot, "/Docs", sep="")
  if(dir.exists(FolderDocs)==FALSE){dir.create(FolderDocs)}
  
  FolderReports =  paste(FolderRoot, "/Reports", sep="")
  if(dir.exists(FolderReports)==FALSE){dir.create(FolderReports)}
  
  FolderDatasets = paste(FolderRoot, "/Datasets", sep="")
  if(dir.exists(FolderScripts)==FALSE){dir.create(FolderScripts)}
  
  FolderCV = paste(FolderDatasets, "/", dataset_name, "/CrossValidation", sep="")
  if(dir.exists(FolderCV)==FALSE){dir.create(FolderCV)}
  
  FolderNamesLabels = paste(FolderDatasets, "/", dataset_name, "/NamesLabels", sep="")
  if(dir.exists(FolderNamesLabels)==FALSE){dir.create(FolderNamesLabels)}
  
  FolderLabelSpace = paste(FolderDatasets, "/", dataset_name, "/LabelSpace", sep="")
  if(dir.exists(FolderLabelSpace)==FALSE){dir.create(LabelSpace)}
  
  FolderTR = paste(FolderCV, "/Tr", sep="")
  if(dir.exists(FolderTR)==FALSE){dir.create(FolderTR)}
  
  FolderTS = paste(FolderCV, "/Ts", sep="")
  if(dir.exists(FolderTS)==FALSE){dir.create(FolderTS)}
  
  FolderVL = paste(FolderCV, "/Vl", sep="")
  if(dir.exists(FolderVL)==FALSE){dir.create(FolderVL)}
  
  FolderSim = paste(FolderRoot, "/Similarities_Results", sep="")
  if(dir.exists(FolderSim)==FALSE){dir.create(FolderSim)}
  
  FolderPlots = paste(FolderRoot, "/Plots", sep="")
  if(dir.exists(FolderPlots)==FALSE){dir.create(FolderPlots)}
  
  retorno$FolderScripts = FolderScripts
  retorno$FolderDocs = FolderDocs
  retorno$FolderReports = FolderReports
  retorno$FolderDatasets = FolderDatasets
  retorno$FolderCV = FolderCV
  retorno$FolderNamesLabels = FolderNamesLabels
  retorno$FolderLabelSpace = FolderLabelSpace
  retorno$FolderTR = FolderTR
  retorno$FolderTS = FolderTS
  retorno$FolderVL = FolderVL
  retorno$FolderSimilarities = FolderSim
  retorno$FolderPlots = FolderPlots
  
  return(retorno)
  gc()
  
}



##################################################################################################