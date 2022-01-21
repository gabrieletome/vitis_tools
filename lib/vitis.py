import sys
import re
import os
import datetime as dt
import time as ti
import zipfile
import lib.filters as filters
import lib.graphic as graphic
import lib.utilities as ut
import glob
from pathlib import Path

autoSaveImg = False
min_frel = 0
list_Genes = []

#Build list of lists of gene
#Return a list of lists
def buildMatrixGenesVitis(listFilter, listFiles):
    global list_Genes
    matrixGenes = []
    listacodici=[]
    extensionFiles = listFiles[0][-4:]
    for f in listFiles:
        if f[-4:] == extensionFiles and (f[-4:] != '.csv' or f[-4:] != '.zip'):
            if f[-4:] == '.zip':
                archive = zipfile.ZipFile(f, 'r')
                fileArchive = archive.namelist()
                for namefilezip in fileArchive:
                    csvText = str(archive.read(namefilezip))
                    if csvText[1] == '"':
                        csvText = csvText.split(r'"')
                    else:
                        csvText = csvText.split('\'')
                    csvText = csvText[1]
                    csvListText = csvText.split(r'\n')
                    csvTemp = open(namefilezip, 'w')
                    for l in csvListText:
                        csvTemp.write(l+'\n')
                    csvTemp.close()
                    listGenes = ut.readFilesVitis(namefilezip)
                    list_Genes = listGenes[1]
                    listGenes = listGenes[0]
                    os.remove(namefilezip)
                    #Filter lists
                    for filter in listFilter:
                        listGenes = applyFilter(listGenes, filter)
                    matrixGenes.append(listGenes)
            elif f[-4:] == '.csv':
                #Read gene files .csv
                listGenes = ut.readFilesVitis(namefilezip)
                list_Genes = listGenes[1]
                listGenes = listGenes[0]
                #Filter lists
                for filter in listFilter:
                    listGenes = applyFilter(listGenes, filter)
                matrixGenes.append(listGenes)
            else:
                #print('cacca')
                #listacodici=[]
                codici=listFiles[0]
                listacodici=codici.split(',')
                file_csv=[]
                ncodici=len(listacodici)
                #print(ncodici)
                #dirname = os.path.dirname(annotated)
                #filename = os.path.join(dirname, 'relative/path/to/file/you/want')
                for x in listacodici:
                    file_csv.append(glob.glob(('../annotated/'+str(x)+'_*')))
                #print(listacodici)
                #print(file_csv)
                #print(listacodici)
                archiviovec=[]
                for x in range(ncodici):
                    if file_csv[x][0] != 0:
                        archiviovec.append(file_csv[x][0])
                #print(archiviovec)
                #////////////////////////////copio il codice da sopra per le zip//////////
                archive = archiviovec
                #print(archive)
                fileArchive = archiviovec
                #print(fileArchive)
                for namefilezip in fileArchive:
                    listGenes = ut.readFilesVitis(namefilezip)
                    list_Genes = listGenes[1]
                    listGenes = listGenes[0]
                    #os.remove(namefilezip)
                    #Filter lists
                    for filter in listFilter:
                        listGenes = applyFilter(listGenes, filter)
                    matrixGenes.append(listGenes)
                    #///////////////////////////////////////////////////////////////////////
        else:
            print('ERROR: FILES HAVE DIFFERENT EXTENSION. File need to have the same extension. All .csv or all .zip')
            sys.exit(-1)
    #return lists of gene lists filtered
    return matrixGenes

#Switch the filter to the correct function
#Return list of genes filtered
def applyFilter(listGenes, filter):
    if filter[0] == '-f' and len(filter) == 2:
        listGenes = filters.filterFrel(listGenes, float(filter[1]))
        global min_frel
        min_frel = float(filter[1])
    elif filter[0] == '-t' and len(filter) >= 2:
        listGenes = filters.filterType(listGenes, filter[1:])
    elif filter[0] == '-a' and len(filter) == 1:
        global autoSaveImg
        autoSaveImg = True
    else:
        print('ERROR: incorrect parameters filter')
        ut.printInfo()
    return listGenes


#Print in output the graphs
def printOutput(coreGraph, graphGenes):
    if not os.path.exists('networkOutput'):
        os.mkdir('networkOutput')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('networkOutput/'+printTime)
    print('Creating directory: \'networkOutput/'+printTime+'\'', flush=True)

    #Pearson correlation
    print('Calculating Pearson Correlation...', flush=True)
    #only on complete graph because all edges in core are in the complete graph
    pearsonComplete = ut.pearsonCorrelation(graphGenes, 'vv_exprdata_2.csv')

    #RETURN CORE
    nameFileCore = 'networkOutput/'+printTime+'/Core_Graph'
    coreGraph = sorted(coreGraph, key=ut.ord)
    ut.printCSV(nameFileCore, coreGraph)

    #RETURN GRAPH
    nameFileCompleteGraph = 'networkOutput/'+printTime+'/Complete_Graph'
    graphGenes = sorted(graphGenes, key=ut.ord)
    ut.printCSV(nameFileCompleteGraph, graphGenes)

    #DRAW GRAPH
    graphic.drawGraph(coreGraph, nameFileCore+'_Circular', pearsonComplete, autoSaveImg, [], 1-min_frel, False)
    graphic.drawGraph(coreGraph, nameFileCore, pearsonComplete, autoSaveImg, [], 1-min_frel, True)
    graphic.drawGraph(graphGenes, nameFileCompleteGraph+'_Circular', pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, False)
    graphic.drawGraph(graphGenes, nameFileCompleteGraph, pearsonComplete, autoSaveImg, list_Genes, 1-min_frel, True)
