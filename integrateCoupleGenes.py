import sys
import re
import os
import datetime as dt
import lib.vitis as vit
import lib.utilities_expansion as utex
import lib.components_graph as comp
import lib.graphic as graphic
import lib.utilities as ut

autoSaveImg = False
completeGraph = False
printDiagram = False
min_frel = 0.1 #default cut frel at 0.1
typeAnalyze = -1 # 0 = frel, 1 = rank, 2 = shared genes
rankCut = -1 #if typeAnalyze == 1; read the number of genes to read
nameDir = 'tmpName'


#Read gene list name files
#Return list of files
def readParameters(input):
    global typeAnalyze
    global rankCut
    listFiles = []
    listFilter = []
    i = 1
    while i < len(input):
        if re.search(r'^-\w$', input[i]):
            filterParam = (input[i],)
            i += 1
            endloop = False
            #loop to take all parameters of a filter
            while i < len(input) and endloop == False:
                #stop if find the next filter or add the params to the tuple of the filter
                if re.search(r'^-', input[i]):
                    endloop = True
                else:
                    filterParam = filterParam + (input[i],)
                    i += 1
            if filterParam[0] == input[i-1] and (filterParam[0] != '-a' and filterParam[0] != '-c' and filterParam[0] != '-e'):
                #return an error if the filter has zero parameters and is not '-a'
                print('ERROR: incorrect number of parameters')
                utex.printInfo()
            else:
                #add tuple to the list
                listFilter.append(filterParam)
        #if find pattern '^-files$', the following parameters are files
        elif re.search(r'^-files$', input[i]):
            i += 1
            while i < len(input):
                listFiles.append(input[i])
                i += 1
        elif re.search(r'^-vitis$', input[i]):
            i += 1
        elif re.search(r'^-rank$', input[i]):
            typeAnalyze = 1
            rankCut = int(input[i+1])
            filterParam = (input[i],input[i+1])
            #add tuple to the list
            listFilter.append(filterParam)
            i += 2
        elif re.search(r'^-frel$', input[i]):
            typeAnalyze = 0
            i += 1
        elif re.search(r'^-shared$', input[i]):
            typeAnalyze = 2
            i += 1
        #if is else, return an error
        else:
            print("ERROR: "+input[i])
            i += 1
    #return an error if the list of file is empty
    if len(listFiles) < 2:
        print('ERROR: no files')
        utex.printInfo()

    #print('FILTER: '+str(listFilter))
    #print('FILES: '+str(listFiles))
    return (listFiles,listFilter)

#print file CSV with edges of graph
def createDir():
    #create directory
    if not os.path.exists('commonGenesOutput'):
        os.mkdir('commonGenesOutput')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('commonGenesOutput/'+printTime)
    global nameDir
    nameDir = 'commonGenesOutput/'+printTime+'/'
    print('Creating directory: \''+nameDir+'\'', flush=True)

#Main function
def main():
    if len(sys.argv) >= 2:
        if sys.argv[1] == '--help':
            utex.printInfo()
        elif sys.argv[1] == '-vitis' and (sys.argv[2] == '-shared' or sys.argv[2] == '-rank' or sys.argv[2] == '-frel'):
            #read parameters
            cmd = readParameters(sys.argv)
            for f in cmd[1]:
                if f[0] == '-f':
                    global min_frel
                    min_frel = float(f[1])
                elif f[0] == '-a':
                    global autoSaveImg
                    autoSaveImg = True
                elif f[0] == '-c':
                    global completeGraph
                    completeGraph = True
                elif f[0] == '-e':
                    global printDiagram
                    printDiagram = True
            #read each file .csv, read first line, check if is our gene
            #If YES -> go further, If NO -> check next file
            #When we have read all expansion of all gene, manage lists
            listCouple = utex.readFiles(cmd[0][0])
            #find common genes in the lists readed
            listCommonGenes = []
            edgesGraph = []
            #read files Vitis
            listFiles = utex.readFilesGenes(cmd[0][1:], listCouple, cmd[1])

            #find common genes
            listCommonGenes = utex.findCommonGenes(listCouple, listFiles)
            if typeAnalyze == 0 or typeAnalyze == 1: #frel or rank
                edgesGraph = utex.buildEdgesFrelRank(listCouple, listFiles)
            elif typeAnalyze == 2: #shared
                edgesGraph = listCommonGenes[0]

            #print CSV with genes share between every gene of LGN
            createDir()
            utex.printCSV(edgesGraph, listCommonGenes[1], nameDir)
            #Draw the Venn diagram, Histogram
            if printDiagram:
                utex.printNumberVenn(listCommonGenes, nameDir)
                graphic.printVenn(listCommonGenes[1], listCouple, nameDir)
                textFiles = utex.readFilesGenes(cmd[0][1:], listCouple, [('-f',min_frel)])
                graphic.printHistogram(edgesGraph, textFiles, nameDir)

            pearsonComplete = []
            for l in edgesGraph:
                #Write genes to download and add the edges
                if completeGraph:
                    #write in a .txt the genes to download
                    listFileToRead = utex.nameAssociateGene(l)
                    if len(listFileToRead) > 0:
                        strToPrint = "To draw edges between associated genes please download lists in file \'"+nameDir+str(edgesGraph.index(l))+"/listToDownload.txt"
                        f = open(nameDir+str(edgesGraph.index(l))+'/listToDownload.txt', 'w')
                        for g in listFileToRead:
                            f.write(str(g) + ', ')
                        f.close()
                        print(strToPrint, flush=True)
                        #wait until they gave the path
                        print('Please write the path of file zip:', flush=True)
                        path = input('>>')
                        print('FILE ZIP: '+path, flush=True)
                        #read file and find edges between them
                        matrixGenes = vit.buildMatrixGenesVitis([('-f', min_frel)], [path])
                        print('Finding edges between genes...')
                        graphGenes = ut.manageDuplicates(ut.buildGraph(matrixGenes))
                        print('Process complete')
                        #save new edges in the list
                        #create dir for each couple of genes
                        nameDirGenes = nameDir+str(edgesGraph.index(l))+'/'
                        #Write file .csv
                        f = open(nameDirGenes+'edges_graph'+'.csv', 'a')
                        f.write('Edges between discovered genes\n')
                        f.write('GeneA,rank,frel,GeneB\n')
                        i = 0
                        while i < len(graphGenes):
                            try:
                                rank = matrixGenes[[u[0] for u in matrixGenes].index(graphGenes[i][0])]
                                rank = rank[[v[1] for v in rank[1:]].index(graphGenes[i][1])]
                            except:
                                rank = matrixGenes[[u[0] for u in matrixGenes].index(graphGenes[i][1])]
                                rank = rank[[v[1] for v in rank[1:]].index(graphGenes[i][0])]
                            #f.write(str(listBioNameUpdate[graphGenes[i][0]])+','+str(rank[0])+','+str(listBioNameUpdate[graphGenes[i][1]])+','+str(graphGenes[i][2])+'\n')
                            f.write(str(graphGenes[i][0])+','+str(rank[0])+','+str(graphGenes[i][2])+','+str(graphGenes[i][1])+'\n')
                            l.append((graphGenes[i][0], rank[0], graphGenes[i][1], graphGenes[i][2]))
                            i += 1
                        f.close()

                #Calculating pearson correlation for each edge
                print('Calculating Pearson correlation...')
                #listForPearson = [((list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(a)]),(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(c)]),d) for (a,b,c,d) in l[1:]]
                #TODEL_GT-001
                # dictConvert = {}
                # for elem in l[1:]:
                #     if elem[0] != 'GT-001':
                #         dictConvert[(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(elem[0])])] = elem[0]
                #     if elem[2] != 'GT-001':
                #         dictConvert[(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(elem[2])])] = elem[2]
                listForPearson = [(a,c,d) for (a,b,c,d) in l[1:]]
                pearsonComplete.append(ut.pearsonCorrelation(listForPearson, 'vv_exprdata_2.csv'))
                #pearson = [(listBioNameUpdate[u],listBioNameUpdate[v],p) for (u,v,p) in tmp]
                #TODEL_GT-001
                #pearson = [(listBioNameUpdate[u],listBioNameUpdate[v],p) for (u,v,p) in tmp]+[((list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(a)]),(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(c)]),1) for (a,b,c,d) in l[1:] if a == 'GT-001' or c == 'GT-001']
                #pearson = [(dictConvert[[w for w in u.split('<BR>') if w in dictConvert.keys()][0]],dictConvert[[w for w in v.split('<BR>') if w in dictConvert.keys()][0]],p) for (u,v,p) in tmp]+[((list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(a)]),(list(listBioNameUpdate.keys())[list(listBioNameUpdate.values()).index(c)]),1) for (a,b,c,d) in l[1:] if a == 'GT-001' or c == 'GT-001']
                #pearsonComplete.append(pearson)
                print('Pearson Correlation done')
            #Draw graph
            graphic.printCommonGraph(edgesGraph, pearsonComplete, 1-min_frel, nameDir, autoSaveImg)
        else:
            print('ERROR: wrong format')
            utex.printInfo()

#Calls the main() function.
if __name__ == '__main__':
    main()
