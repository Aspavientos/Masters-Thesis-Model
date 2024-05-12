# Information ----
# Introducing new datasets from GEO
# Author: Diego Rodríguez Esperante
# Date of creation: 24/03/2024
# Last edited: 24/03/2024

# Packages ----
require(RCy3)
require(rstudioapi)
require(utils)
require(stringr)
require(beepr)

# Connect to Cytoscape ----
cytoscapePing()
cytoscapeVersionInfo()

closeSession(FALSE)

# Set path ----
workingdir = dirname(rstudioapi::getSourceEditorContext()$path)

# Select model to open ----
setwd(workingdir)
network_file = choose.files(getwd(), caption = 'Choose a network file')
importNetworkFromFile(network_file)

nodeData_file = paste0(dirname(network_file), '/',
                       tools::file_path_sans_ext(basename(network_file)), '_NodeData.csv')
nodeData = read.csv(nodeData_file)

loadTableData(nodeData, data.key.column = 'Node')
deleteTableColumn('Node')

# Apply styles ----
setwd(workingdir)
styles_folder = choose.dir(getwd(), caption = 'Choose the folder with style information')

styles_ext = tools::file_ext(dir(styles_folder))
for(i in 1:length(styles_ext)){
  if(styles_ext[i] == 'xml'){
    importVisualStyles(paste0(styles_folder, '\\', dir(styles_folder)[i]))
  }else if(styles_ext[i] == 'txt'){
    nodeLocs = read.table(paste0(styles_folder, '\\', dir(styles_folder)[i]), sep = '\t', header = TRUE)
    loadTableData(nodeLocs, data.key.column = 'Node')
    deleteTableColumn('Node')
  }
}

setVisualStyle('COBRA Expression')

nodeTable = getTableColumns()
hiddennodes = nodeTable$name[is.na(nodeTable$x)]
if(length(hiddennodes)>0){
  hideNodes(nodeTable$name[is.na(nodeTable$x)])
}

hideAllPanels()
fitContent()

# Select data to represent ----
setwd(workingdir)
FBA_files = choose.files(getwd(), paste0('Choose the FBA files (', tools::file_path_sans_ext(basename(network_file)), ')'))

for(i in 1:length(FBA_files)){
  FBAData = read.csv(FBA_files[i])
  loadTableData(FBAData, data.key.column = 'Reactions')
  deleteTableColumn('Reactions')
}

# deleteAnnotation(sapply(getAnnotationList(), '[[', 'uuid'))

addAnnotationText(text = 'Cohort', name = 'Cohort',
                  x.pos = -1700, y.pos = 250, fontSize = 50, fontFamily = 'Agency FB', fontStyle = 'bold')

addAnnotationText(text = 'Data', name = 'Data',
                  x.pos = -1700, y.pos = 300, fontSize = 40, fontFamily = 'Agency FB')


text_uuid = c(0, 0)
text_uuid[1] = getAnnotationList()[[2]][['uuid']]
text_uuid[2] = getAnnotationList()[[1]][['uuid']]

# Cycle through values ----
nodeTable = getTableColumns()
FBA_cols = (which(colnames(nodeTable) == 'y')+1):ncol(nodeTable)

maxmins_cols = as.data.frame(str_split(colnames(nodeTable)[FBA_cols], '\\.', simplify = TRUE), stringsAsFactors = TRUE)
maxmins_cols = cbind(maxmins_cols, data.frame(max = matrix(0, nrow(maxmins_cols), 1), min = matrix(0, nrow(maxmins_cols), 1)))
colnames(maxmins_cols)[1:2] = c('Data', 'Cohort')

for(i in 1:nrow(maxmins_cols)){
  if(maxmins_cols$Data[i] == 'UB'){
    maxmins_cols$max[i] = unique(sort(nodeTable[,FBA_cols[i]], decreasing = TRUE))[2]
    maxmins_cols$min[i] = -1e4
  }else{
    maxmins_cols$max[i] = max(nodeTable[,FBA_cols[i]], na.rm = TRUE)
    maxmins_cols$min[i] = min(nodeTable[,FBA_cols[i]], na.rm = TRUE)
  }
  
}

# mapVisualProperty
setwd(workingdir)
images_folder = choose.dir(caption = 'Select folder to save images')

numim = 1
for(i in FBA_cols){
  message(paste('Exporting', as.character(numim), 'of', length(FBA_cols), 'images...'))
  
  splitname = str_split(colnames(nodeTable)[i], '\\.')[[1]]
  updateAnnotationText(splitname[2], text_uuid[1],
                       x.pos = -1700, y.pos = 250, fontSize = 50, fontFamily = 'Agency FB', fontStyle = 'bold')
  updateAnnotationText(splitname[1], text_uuid[2],
                       x.pos = -1700, y.pos = 300, fontSize = 40, fontFamily = 'Agency FB')
  if(splitname[1] == 'Overr'){
    visStyle = 'COBRA Overrepresentation'
    colors = c('#FFFFFF', '#6BAED6', '#08306B')
    maxmins = c(min(nodeTable[i], na.rm = TRUE),
                max(nodeTable[i], na.rm = TRUE)/2,
                max(nodeTable[i], na.rm = TRUE))
  }else{
    visStyle = 'COBRA Expression'
    colors = c('#00C16A', '#FFFFFF', '#5800BC')
    maxmins = c(min(maxmins_cols$min[maxmins_cols$Data == splitname[1]]),
                0,
                max(maxmins_cols$max[maxmins_cols$Data == splitname[1]]))
  }
  
  setVisualStyle(visStyle)
  
  visMapping = mapVisualProperty('Node Border Paint', colnames(nodeTable)[i], 'c', maxmins, colors)
  updateStyleMapping(visStyle, visMapping)
  
  subfolder_path = paste0(images_folder, '\\', splitname[1], '\\')
  
  if(!dir.exists(subfolder_path)){
    dir.create(subfolder_path)
  }
  
  exportPNG(paste0(subfolder_path, splitname[1], '_', splitname[2]),
            transparentBackground = TRUE, zoom = 500)
  
  numim = numim + 1
}

message('Done!')
beep()
