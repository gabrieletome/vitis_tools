import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib as mpl
import networkx as nx
import numpy as np
from operator import itemgetter
from matplotlib_venn import venn3, venn2
import lib.components_graph as comp
import lib.utilities_expansion as utex
import lib.utilities as ut
import lib.venn as vennD
import re
import os
import random
from networkx.readwrite import json_graph
import json
import sys

#trovo l indice di un nodo in listlinename
def nodefounder(a,b,c):
            quantograndelistlinename=len(a)
            for y in range(quantograndelistlinename):
                nododatrovare=b[c]
                indicetrovato=y
                #global indicenodotrovato
                if (nododatrovare==a[y][0]):
                    #indicetrovato=y
                    print(indicetrovato)
                    return indicetrovato

#Draw general graph
def drawGraph(net, namefile, pearson, autoSaveImg, list_Genes, range_frel, typePrint):
    #range to divide type of edges
    step1_range = round((range_frel/3)+1-range_frel, 2)
    step2_range = round((range_frel/3*2)+1-range_frel, 2)
    #set name of window and change parameters of save image
    namefiles = namefile.split('/')
    plt.figure(namefiles[2])
    mpl.rcParams['savefig.dpi'] = 700
    mpl.rcParams['savefig.directory'] = namefiles[0]+'/'+namefiles[1]
    #convert name genes to number
    idNode = {}
    i = 1
    edges = []
    newPearson = []
    for k in net:
        if k[0] not in idNode:
            idNode[k[0]] = i
            i += 1
        if k[1] not in idNode:
            idNode[k[1]] = i
            i += 1
    for k in net:
        edges.append((idNode[k[0]], idNode[k[1]], k[2]))
    for k in pearson:
        try:
            newPearson.append((idNode[k[0]], idNode[k[1]], k[2]))
        except:
            pass #edge of complete graph not in core

    for k in list_Genes:
        if k not in idNode.keys():
            edges.append((i, i, 1))
            newPearson.append((i, i, 1))
            idNode[k] = i
            i += 1

    #order id based on name gene
    oldidNode = idNode.copy()
    tmpKeys = sorted(idNode.keys())
    for k in tmpKeys:
        idNode[k] = tmpKeys.index(k)+1
    i = 0
    while i < len(edges):
        tmp = edges[i]
        edges[i] = (idNode[list(oldidNode.keys())[list(oldidNode.values()).index(tmp[0])]], idNode[list(oldidNode.keys())[list(oldidNode.values()).index(tmp[1])]], tmp[2])
        i += 1
    i = 0
    while i < len(newPearson):
        tmp = newPearson[i]
        newPearson[i] = (idNode[list(oldidNode.keys())[list(oldidNode.values()).index(tmp[0])]], idNode[list(oldidNode.keys())[list(oldidNode.values()).index(tmp[1])]], tmp[2])
        i += 1

    G = nx.Graph()
    for e in edges:
        G.add_edge(e[0], e[1], weight=e[2])





    #print the connected components of the complete graph
    if len(list_Genes) > 0 and not typePrint:
        subgraphs = sorted(list(G.subgraph(c) for c in nx.connected_components(G)), key=len, reverse=True)
        f = open(namefiles[0]+"/"+namefiles[1]+'/complete_graph_connected_components.txt', 'w')
        i = 0
        while i < len(subgraphs):
            s = subgraphs[i]
            if len(s.nodes()) > 1:
                f.write('Node subgraph '+str(i+1)+':\n')
                lNodes = []
                for n in s.nodes():
                    lNodes.append((list(idNode.keys())[list(idNode.values()).index(n)]))
                f.write(str(lNodes)+'\n')
            i += 1
        f.close()
        #print('Create: '+namefiles[0]+"/"+namefiles[1]+'/complete_graph_connected_components.txt', flush=True)

    ePos = [(u, v) for (u, v, p) in newPearson if p >= 0]
    eNeg = [(u, v) for (u, v, p) in newPearson if p < 0]

    estrongPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > step2_range and ((u, v) in ePos or (v, u) in ePos)]
    estrongNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > step2_range and ((u, v) in eNeg or (v, u) in eNeg)]
    emediumPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step2_range and d['weight'] > step1_range and ((u, v) in ePos or (v, u) in ePos)]
    emediumNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step2_range and d['weight'] > step1_range and ((u, v) in eNeg or (v, u) in eNeg)]
    eweakPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step1_range and ((u, v) in ePos or (v, u) in ePos)]
    eweakNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step1_range and ((u, v) in eNeg or (v, u) in eNeg)]
    # positions for all nodes
    pos = comp.layout_many_components(G, component_layout_func=nx.layout.circular_layout)
    if typePrint:
        pos = comp.layout_many_components(G, component_layout_func=nx.layout.spring_layout)
        #print Degree of graph
        printDegree(G, idNode, namefile)
    # nodes
    n_size = 150
    f_size = 8
    e_size = 2
    if len(idNode.keys()) > 100:
        n_size = 7
        f_size = 1
        e_size = 0.2
    elif len(idNode.keys()) >= 50 and len(idNode.keys()) < 100:
        n_size = 15
        f_size = 2
        e_size = 0.5
    elif len(idNode.keys()) >= 30 and len(idNode.keys()) < 50:
        n_size = 50
        f_size = 4
        e_size = 1
    elif len(idNode.keys()) >= 15 and len(idNode.keys()) < 30:
        n_size = 150
        f_size = 8
        e_size = 1.5
    #draw node vitis
    nx.draw_networkx_nodes(G, pos, node_size=n_size, node_color='#C7E7EB')
    # edges
    nx.draw_networkx_edges(G, pos, edgelist=estrongPos, width=e_size, edge_color='black')
    nx.draw_networkx_edges(G, pos, edgelist=emediumPos, width=e_size/2, edge_color='black', style='dashed')
    nx.draw_networkx_edges(G, pos, edgelist=eweakPos, width=e_size/2, edge_color='black', style='dotted')
    nx.draw_networkx_edges(G, pos, edgelist=estrongNeg, width=e_size, edge_color='red')
    nx.draw_networkx_edges(G, pos, edgelist=emediumNeg, width=e_size/2, edge_color='red', style='dashed')
    nx.draw_networkx_edges(G, pos, edgelist=eweakNeg, width=e_size/2, edge_color='red', style='dotted')
    # labels
    nx.draw_networkx_labels(G, pos, font_size=f_size, font_family='sans-serif')
    plt.axis('off')
    #legend
    black_line = mlines.Line2D([], [], linewidth=2, color='black', label='frel > '+str(step2_range))
    blue_line = mlines.Line2D([], [], linewidth=1, color='black', label=str(step1_range)+' < frel <= '+str(step2_range), linestyle='dashed')
    red_line = mlines.Line2D([], [], linewidth=1, color='black', label='frel <= '+str(step1_range), linestyle='dotted')
    pos_line = mlines.Line2D([], [], linewidth=1, color='black', label='Pearson corr. >= 0')
    neg_line = mlines.Line2D([], [], linewidth=1, color='red', label='Pearson corr. < 0')
    if len(eNeg) > 0:
        textLegend = [black_line, blue_line, red_line, pos_line, neg_line]
    else:
        textLegend = [black_line, blue_line, red_line]
    #Write legend ID-->GENE
    nameGenes = idNode.keys()
    dictStrToWrite = {}
    #Read information of Vitis genes
    f = open('import_doc/AnnVitisNetV1_V3.csv', 'r')
    text = f.readlines()
    listLineName = []
    i = 1
    while i < len(text):
        listLineName.append(text[i].split(','))
        i += 1
    for k in nameGenes:
        if k in [u[0].upper() for u in listLineName]:
            index = [u[0].upper() for u in listLineName].index(k)
            u = listLineName[index]
            dictStrToWrite[k] = str(u[0])+','+str(u[1])+','+str(u[2])+','+str(u[3])+','+str(u[4])+','+str(u[5])+','+str(u[6])
        else:
            dictStrToWrite[k] = str(k)
    fileOut = namefile.split("Graph")[0]+'graph_legend_ID_NAME.csv'
    #print('LEGEND IN: \''+fileOut+'\'')
    f = open(fileOut, 'w')
    f.write('ID in graph,'+text[0].split(',')[0]+','+text[0].split(',')[1]+','+text[0].split(',')[2]+','+text[0].split(',')[3]+','+text[0].split(',')[4]+','+text[0].split(',')[5]+','+text[0].split(',')[6]+',\n')
    for k in nameGenes:
        if dictStrToWrite[k].endswith('\n'):
            f.write(str(idNode[k])+','+dictStrToWrite[k])
        else:
            f.write(str(idNode[k])+','+dictStrToWrite[k]+'\n')
    f.close()

    plt.legend(handles=textLegend, fontsize = 'xx-small').set_draggable(True)
    #autoSave PNG or show
    if autoSaveImg:
        plt.savefig(namefile+'.png')
        #print('Create: \''+namefile+'.png\'', flush=True)
    else:
        plt.show()





    id_nodi_acasissimo=list(G.nodes)

    #Create json file of graph
    if typePrint:
        H = nx.relabel_nodes(G, dict((v,k) for k,v in idNode.items()))
        #print('Create: \''+namefile+'.json\'', flush=True)
        # with open(namefile+'.json', 'w', encoding='utf-8') as f:
        #     json.dump(json_graph.node_link_data(G), f, ensure_ascii=False, indent=4)
        #     f.close



        #JSON COMPATIBILE CON CYTOSCAPE//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        #H = G
        #print('prontone')
        #H=nx.relabel_nodes(H,dict)
        lista_deinodi=list(H.nodes)
        #print(lista_deinodi)
        lista_degliarchi=list(G.edges)
        #print(lista_degliarchi)
        n_Archi=len(lista_degliarchi)
        n_nodi=len(lista_deinodi)
        lista_nodi_Archi=['nodes','edges']
        id_nodes=[]
        id_edges=[]


        #///liste degli id/////////////////////////////
        for x in range(n_nodi):
            id_nodes.append(x)
            id_nodes[x]=x
        #print(id_nodes)
        for x in range(n_Archi):
            id_edges.append(x)
            id_edges[x]=x
        #print(id_edges)
        listof_id=[]
        listof_id.append(id_nodes)
        listof_id.append(id_edges)
        #print(listof_id)



        #//creo l elemento data con gli attributi(dizionario di dizionari)//////////////////////////////////

        data_nodes=[]
        data_edges=[]
        dizionario_nested_nodi={}
        dizionario_nested_archi={}

        for x in range(n_nodi):
            #print(lista_deinodi[x])
            i1={"id":lista_deinodi[x]}                                    #qui posso inserire attributi nuovi per i nodi
            i2={"label":lista_deinodi[x]}
            i3={"vit3":'inserire qui vit3'}
            i4={"vit2":"inserire qui vit2"}

            dizionario_nested_nodi={"data":(i1,i2,i3,i4)}


            data_nodes.append(x)
            data_nodes[x]=dizionario_nested_nodi
        #print(data_nodes)
        n_attributi_nodi=4
        n_attributi_archi=4
        for x in range(n_Archi):
            #print(lista_degliarchi[x])
            i1={"id":x}                                                 #qui posso inserire attributi nuovi per gli archi
            i2={"source":lista_degliarchi[x][0]}
            i3={"target":lista_degliarchi[x][1]}
            i4={"colourId":'boh'}
            dizionario_nested_archi={"data":[i1,i2,i3,i4]}


            data_edges.append(x)
            data_edges[x]=dizionario_nested_archi
        #print(data_edges)



        #//////dizionario tra id e data//////////
        id_data_n_diz={}
        #id_data_n_diz=dict(zip(id_nodes,data_nodes))
        #print(id_data_n_diz)

        id_data_e_diz={}
        #id_data_e_diz=dict(zip(id_edges,data_edges))
        #print(id_data_e_diz)

        sumofcosette=[]
        sumofcosette.append(data_nodes)
        sumofcosette.append(data_edges)

        #Write legend ID-->GENE

        #print(listLineName)

        #///////aggiungo id-data a nodes e edges//////////////////////////////////////////
        dizionariobello={}
        dizionariobello=dict(zip(lista_nodi_Archi,sumofcosette))

        spaghetti=list(G.edges)
        quantograndelistlinename=len(listLineName)
        dizzi={}
        listabela=[]
        for k in nameGenes:
            if k in [u[0].upper() for u in listLineName]:
                index = [u[0].upper() for u in listLineName].index(k)
                u = listLineName[index]
                listabela=[u[0],u[1],u[2],u[3],u[4],u[5],u[6]]
                dizzi[k] = listabela
            else:
                dizzi[k] = str(k)+"\n"

        #/////creo il file txt per poterlo modificere come testo e poi passarlo in un json///
        open(namefile+'temp.json', 'w', encoding='utf-8')
        with open(namefile+'temp.json', 'a', encoding='utf-8') as f:
            # f.write("lol"+dizionariobello['nodes'][1]['data'][1]['label'])
            #f.write('{"elements":')
            #spaghetti=list(G.edges)
            f.write('{"elements":')
            f.write('{')
            f.write('"nodes": [')
            #indicenodotrovato
            for x in range(n_nodi):
                nodoattuale=lista_deinodi[x]
                f.write('{')
                f.write('"data":')
                f.write('  {')
                f.write('"id":"'+str(id_nodi_acasissimo[x])+'",')
                f.write('"label":"'+ lista_deinodi[x]+'",')
                f.write('"functann":"'+dizzi[nodoattuale][1]+'",')
                f.write('"nameEC":"'+dizzi[nodoattuale][2]+'",')
                f.write('"gn":"'+dizzi[nodoattuale][3]+'",')
                f.write('"vit3":"'+dizzi[nodoattuale][4]+'",')
                f.write('"nw1":"'+dizzi[nodoattuale][5]+'",')
                f.write('"nw2":"'+dizzi[nodoattuale][6]+'",')
                f.write('"nodec":"0"')
                f.write('}')
                f.write('},')
            f.write('],')
            f.write('"edges": [')
            for x in range(n_Archi):
                f.write('{')
                f.write('"data":')
                f.write('  {')
                f.write('"id":"e'+ str(x)+'",')
                f.write('"source":"'+str(lista_degliarchi[x][0])+'",')
                f.write('"target":"'+str(lista_degliarchi[x][1])+'",')
                f.write('"weight":"'+'",')
                f.write('"codec":')

                if (spaghetti[x][0],spaghetti[x][1]) in estrongPos:
                    f.write('"1"')
                elif (spaghetti[x][0],spaghetti[x][1]) in estrongNeg:
                    f.write('"2"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumPos:
                    f.write('"3"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumNeg:
                    f.write('"4"')
                elif  (spaghetti[x][0],spaghetti[x][1])in eweakPos:
                    f.write('"5"')
                else:
                    f.write('"6"')

                f.write('}')
                f.write('},')
            f.write(']')
            f.write('}')
            f.write('}')
            f.close

            #da printare su terminale//////
            '''
            spaghetti=list(G.edges)
            print('{')
            print('\t'+'"nodes": [')
            for x in range(n_nodi):
                print('\t\t\t'+'   {')
                print('\t\t\t\t'+'"data":')
                print('\t\t\t\t'+'  {')
                print('\t\t\t\t\t'+'"id":"'+ str(x)+'",')
                print('\t\t\t\t\t'+'"label":"'+ lista_deinodi[x]+'",')
                print('\t\t\t\t\t'+'"vit3":"'+'",')
                print('\t\t\t\t\t'+'"vit2":"'+'"')
                print('\t\t\t\t'+'}')
                print('\t\t\t     },')
            print('\t\t\t],')
            print('\t'+'"edges": [')
            for x in range(n_Archi):
                print('\t\t\t'+'   {')
                print('\t\t\t\t'+'"data":')
                print('\t\t\t\t'+'  {')
                print('\t\t\t\t\t'+'"id":"'+ str(x)+'",')
                print('\t\t\t\t\t'+'"source":"'+lista_degliarchi[x][0]+'",')
                print('\t\t\t\t\t'+'"target":"'+lista_degliarchi[x][1]+'",')
                print('\t\t\t\t\t'+'"wheight":"'+'",')
                print('\t\t\t\t\t'+'"colourId":')
                if (spaghetti[x][0],spaghetti[x][1]) in estrongPos:
                    print('"estrongpos"')
                elif (spaghetti[x][0],spaghetti[x][1]) in estrongNeg:
                    print('"estrongneg"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumPos:
                    print('"emediumpos"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumNeg:
                    print('"emediumneg"')
                elif  (spaghetti[x][0],spaghetti[x][1])in eweakPos:
                    print('"eweakpos"')
                else:
                    print('"eweakneg"')
                print('\t\t\t\t'+'}')
                print('\t\t\t     },')
            print('\t\t\t]')
            print('}')
            print(spaghetti)
            #print(lista_deinodi)
            #print(lista_degliarchi)
            print(estrongNeg)
            #print(ePos)
            #print(eNeg)
            '''

        #modifico json///////////////////////////////////
        polenta=[]
        with open(namefile+'temp.json', 'r', encoding='utf-8') as f:
            polenta=f.read()
            #print(polenta)
        polentaold=polenta
        pattern1= re.compile(r'\,\]')
        pattern2= re.compile(r'\,')
        pattern3=re.compile(r'\[')
        pattern4=re.compile(r'\]')
        pattern5=re.compile(r'\{')
        pattern6=re.compile(r'\}')
        pattern7=re.compile(r'\n\"')

        if pattern1.search(polenta):
            polenta=pattern1.sub(r']',polenta)
        if pattern7.search(polenta):
            polenta=pattern7.sub(r'"',polenta)
        if pattern2.search(polenta):
            polenta=pattern2.sub(r',\n',polenta)
        if pattern3.search(polenta):
            polenta=pattern3.sub(r'[\n',polenta)
        if pattern4.search(polenta):
            polenta=pattern4.sub(r']\n',polenta)
        if pattern5.search(polenta):
            polenta=pattern5.sub(r'{\n',polenta)
        if pattern6.search(polenta):
            polenta=pattern6.sub(r'}\n',polenta)


        stampare=polenta
        f=open(namefile+'final.json', 'w', encoding='utf-8')
        f.write(stampare)
        f.close
        #print(polenta)

        '''
        filetesto = open(namefile+"yo.txt","r")
        polentaold=filetesto.read()
        polenta=polentaold
        pattern= re.compile(r',]')
        with open("yo.txt") as f:
            for line in f:
                if pattern.search(line):
                    polenta=pattern.sub(r']',polenta)
        '''
        os.remove(namefile+'temp.json')
    #/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #Clean graph and pyplot
    G.clear()
    plt.clf()
    plt.close()



#Draw graph of expansion LGN
def printCommonGraph(listCommonGenes, pearsonComplete, range_frel, nameDir, autoSaveImg):
    #range to divide type of edges
    step1_range = round((range_frel/3)+1-range_frel, 2)
    step2_range = round((range_frel/3*2)+1-range_frel, 2)

    for l in listCommonGenes:
        #nameFile
        nameF = utex.buildNamefile(l)
        namefile = nameDir+str(listCommonGenes.index(l))+'/graphGenes'
        #set name of window and change parameters of save image
        plt.figure('graph')
        mpl.rcParams['savefig.directory'] = nameDir+nameF
        mpl.rcParams['savefig.dpi'] = 700
        #build edges of graph
        net = utex.buildEdges(l)
        pearson = pearsonComplete[listCommonGenes.index(l)]
        #convert name genes to number
        idNode = {}
        i = 1
        edges = []
        newPearson = []
        for k in net:
            if k[0] not in idNode:
                idNode[k[0]] = i
                i += 1
            if k[1] not in idNode:
                idNode[k[1]] = i
                i += 1
        for k in net:
            edges.append((idNode[k[0]], idNode[k[1]], k[2]))
        for k in pearson:
            try:
                newPearson.append((idNode[k[0]], idNode[k[1]], k[2]))
            except:
                pass #edge of complete graph not in core

        #Set color of nodes. Genes of LGN -> RED, associated genes -> BLUE
        dict_isoColor = {}
        isoform_done = {}
        for g in edges:
            #print((list(idNode.keys())[list(idNode.values()).index(g[0])]))
            if (list(idNode.keys())[list(idNode.values()).index(g[0])]) in l[0] and (list(idNode.keys())[list(idNode.values()).index(g[0])]) not in isoform_done.keys():
                dict_isoColor[g[0]] = "#FF3232"
                isoform_done[(list(idNode.keys())[list(idNode.values()).index(g[0])])] = g[0]
            elif (list(idNode.keys())[list(idNode.values()).index(g[0])]) not in isoform_done.keys():
                dict_isoColor[g[0]] = '#48B3FF'
                isoform_done[(list(idNode.keys())[list(idNode.values()).index(g[0])])] = g[0]
            if (list(idNode.keys())[list(idNode.values()).index(g[1])]) in l[0] and (list(idNode.keys())[list(idNode.values()).index(g[1])]) not in isoform_done.keys():
                dict_isoColor[g[1]] = "#FF3232"
                isoform_done[(list(idNode.keys())[list(idNode.values()).index(g[1])])] = g[1]
            elif (list(idNode.keys())[list(idNode.values()).index(g[1])]) not in isoform_done.keys():
                dict_isoColor[g[1]] = '#48B3FF'
                isoform_done[(list(idNode.keys())[list(idNode.values()).index(g[1])])] = g[1]

        #add edges to Graph
        G = nx.Graph()
        for e in edges:
            G.add_edge(e[0], e[1], weight=e[2])

        #divide edges by weight
        ePos = [(u, v) for (u, v, p) in newPearson if p >= 0]
        eNeg = [(u, v) for (u, v, p) in newPearson if p < 0]
        estrongPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > step2_range and ((u, v) in ePos or (v, u) in ePos)]
        estrongNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] > step2_range and ((u, v) in eNeg or (v, u) in eNeg)]
        emediumPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step2_range and d['weight'] > step1_range and ((u, v) in ePos or (v, u) in ePos)]
        emediumNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step2_range and d['weight'] > step1_range and ((u, v) in eNeg or (v, u) in eNeg)]
        eweakPos = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step1_range and ((u, v) in ePos or (v, u) in ePos)]
        eweakNeg = [(u, v) for (u, v, d) in G.edges(data=True) if d['weight'] <= step1_range and ((u, v) in eNeg or (v, u) in eNeg)]
        pos = comp.layout_many_components(G, component_layout_func=nx.layout.spring_layout)

        # size of nodes, text and edges
        n_size = 150
        f_size = 8
        e_size = 2
        if len(idNode.keys()) > 100:
            n_size = 7
            f_size = 1
            e_size = 0.2
        elif len(idNode.keys()) >= 50 and len(idNode.keys()) < 100:
            n_size = 15
            f_size = 2
            e_size = 0.5
        elif len(idNode.keys()) >= 30 and len(idNode.keys()) < 50:
            n_size = 50
            f_size = 4
            e_size = 1
        elif len(idNode.keys()) >= 15 and len(idNode.keys()) < 30:
            n_size = 100
            f_size = 6
            e_size = 1.5

        #draw nodes
        vector_isoColor = []
        if len(dict_isoColor) > 0:
            for node in G.nodes():
                vector_isoColor.append(dict_isoColor[node])
        nx.draw_networkx_nodes(G, pos, node_size=n_size, node_color=vector_isoColor)
        vector_isoColor = []
        if len(dict_isoColor) > 0:
            for node in G.nodes():
                vector_isoColor.append(dict_isoColor[node])
        nx.draw_networkx_nodes(G, pos, node_size=n_size, node_color=vector_isoColor)
        # edges
        nx.draw_networkx_edges(G, pos, edgelist=estrongPos, width=e_size, edge_color='black')
        nx.draw_networkx_edges(G, pos, edgelist=emediumPos, width=e_size/2, edge_color='black', style='dashed')
        nx.draw_networkx_edges(G, pos, edgelist=eweakPos, width=e_size/2, edge_color='black', style='dotted')
        nx.draw_networkx_edges(G, pos, edgelist=estrongNeg, width=e_size, edge_color='red')
        nx.draw_networkx_edges(G, pos, edgelist=emediumNeg, width=e_size/2, edge_color='red', style='dashed')
        nx.draw_networkx_edges(G, pos, edgelist=eweakNeg, width=e_size/2, edge_color='red', style='dotted')
        # labels
        nx.draw_networkx_labels(G, pos, font_size=f_size, font_family='sans-serif')

        plt.axis('off')
        #plt.title("Graph LGN with related genes")
        #legend
        black_line = mlines.Line2D([], [], linewidth=2, color='black', label='frel > '+str(step2_range))
        blue_line = mlines.Line2D([], [], linewidth=1, color='black', label=str(step1_range)+' < frel <= '+str(step2_range), linestyle='dashed')
        red_line = mlines.Line2D([], [], linewidth=1, color='black', label='frel <= '+str(step1_range), linestyle='dotted')
        pos_line = mlines.Line2D([], [], linewidth=1, color='black', label='Pearson corr. >= 0')
        neg_line = mlines.Line2D([], [], linewidth=1, color='red', label='Pearson corr. < 0')
        if len(eNeg) > 0:
            textLegend = [black_line, blue_line, red_line, pos_line, neg_line]
        else:
            textLegend = [black_line, blue_line, red_line]
        circleLGN = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor="#FF3232", label='Genes in LGN')
        circleNewNodes = mlines.Line2D(range(1), range(1), color="white", marker='o', markerfacecolor='#48B3FF', label='Discovered genes')
        textLegend.append(circleLGN)
        textLegend.append(circleNewNodes)
        #Write legend ID-->GENE
        nameGenes = idNode.keys()
        dictStrToWrite = {}
        #Read information of Vitis genes
        f = open('import_doc/AnnVitisNetV1_V3.csv', 'r')
        text = f.readlines()
        listLineName = []
        i = 1
        while i < len(text):
            listLineName.append(text[i].split(','))
            i += 1
        for k in nameGenes:
            if k in [u[0].upper() for u in listLineName]:
                index = [u[0].upper() for u in listLineName].index(k)
                u = listLineName[index]
                dictStrToWrite[k] = str(u[0])+','+str(u[1])+','+str(u[2])+','+str(u[3])+','+str(u[4])+','+str(u[5])+','+str(u[6])
            else:
                dictStrToWrite[k] = str(k)+"\n"
        fileOut = namefile.split("graph")[0]+'graph_legend_ID_NAME.csv'
        print('LEGEND IN: \''+fileOut+'\'')
        f = open(fileOut, 'w')
        f.write('ID in graph,'+text[0].split(',')[0]+','+text[0].split(',')[1]+','+text[0].split(',')[2]+','+text[0].split(',')[3]+','+text[0].split(',')[4]+','+text[0].split(',')[5]+','+text[0].split(',')[6])
        for k in nameGenes:
            if dictStrToWrite[k].endswith('\n'):
                f.write(str(idNode[k])+','+dictStrToWrite[k])
            else:
                f.write(str(idNode[k])+','+dictStrToWrite[k]+'\n')
        f.close()

        #textLegend.append(mlines.Line2D([], [], label='in \''+(fileOut.split('/'))[-1]+'\'', visible=False))
        plt.legend(handles=textLegend, fontsize = 'xx-small').set_draggable(True)
        #autoSave PNG or show
        if autoSaveImg:
            namefile = namefile.replace("<", "_")
            namefile = namefile.replace(">", "_")
            plt.savefig(namefile+'.png')
            print('Create: \'graphGenes.png\'', flush=True)
        else:
            plt.show()

        #print Degree of graph
        printDegree(G, idNode, namefile)



        #Create json file of graph

        H = nx.relabel_nodes(G, dict((v,k) for k,v in idNode.items()))
        # #print('Create: \''+namefile+'.json\'', flush=True)
        # with open(namefile+'.json', 'w', encoding='utf-8') as f:
        #     json.dump(json_graph.node_link_data(G), f, ensure_ascii=False, indent=4)
        #     f.close


        id_nodi_acasissimo=list(G.nodes)
        #JSON COMPATIBILE CON CYTOSCAPE//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        #H = G
        #print('prontone')
        #H=nx.relabel_nodes(H,dict)
        lista_deinodi=list(H.nodes)
        #print(lista_deinodi)
        lista_degliarchi=list(G.edges)
        #print(lista_degliarchi)
        n_Archi=len(lista_degliarchi)
        n_nodi=len(lista_deinodi)
        lista_nodi_Archi=['nodes','edges']
        id_nodes=[]
        id_edges=[]


        #///liste degli id/////////////////////////////
        for x in range(n_nodi):
            id_nodes.append(x)
            id_nodes[x]=x
        #print(id_nodes)
        for x in range(n_Archi):
            id_edges.append(x)
            id_edges[x]=x
        #print(id_edges)
        listof_id=[]
        listof_id.append(id_nodes)
        listof_id.append(id_edges)
        #print(listof_id)



        #//creo l elemento data con gli attributi(dizionario di dizionari)//////////////////////////////////

        data_nodes=[]
        data_edges=[]
        dizionario_nested_nodi={}
        dizionario_nested_archi={}

        for x in range(n_nodi):
            #print(lista_deinodi[x])
            i1={"id":lista_deinodi[x]}                                    #qui posso inserire attributi nuovi per i nodi
            i2={"label":lista_deinodi[x]}
            i3={"vit3":'inserire qui vit3'}
            i4={"vit2":"inserire qui vit2"}

            dizionario_nested_nodi={"data":(i1,i2,i3,i4)}


            data_nodes.append(x)
            data_nodes[x]=dizionario_nested_nodi
        #print(data_nodes)
        n_attributi_nodi=4
        n_attributi_archi=4
        for x in range(n_Archi):
            #print(lista_degliarchi[x])
            i1={"id":x}                                                 #qui posso inserire attributi nuovi per gli archi
            i2={"source":lista_degliarchi[x][0]}
            i3={"target":lista_degliarchi[x][1]}
            i4={"colourId":'boh'}
            dizionario_nested_archi={"data":[i1,i2,i3,i4]}


            data_edges.append(x)
            data_edges[x]=dizionario_nested_archi
        #print(data_edges)



        #//////dizionario tra id e data//////////
        id_data_n_diz={}
        #id_data_n_diz=dict(zip(id_nodes,data_nodes))
        #print(id_data_n_diz)

        id_data_e_diz={}
        #id_data_e_diz=dict(zip(id_edges,data_edges))
        #print(id_data_e_diz)

        sumofcosette=[]
        sumofcosette.append(data_nodes)
        sumofcosette.append(data_edges)

        #Write legend ID-->GENE

        #print(listLineName)

        #///////aggiungo id-data a nodes e edges//////////////////////////////////////////
        dizionariobello={}
        dizionariobello=dict(zip(lista_nodi_Archi,sumofcosette))

        spaghetti=list(G.edges)
        quantograndelistlinename=len(listLineName)
        dizzi={}
        listabela=[]
        for k in nameGenes:
            if k in [u[0].upper() for u in listLineName]:
                index = [u[0].upper() for u in listLineName].index(k)
                u = listLineName[index]
                listabela=[u[0],u[1],u[2],u[3],u[4],u[5],u[6]]
                dizzi[k] = listabela
            else:
                dizzi[k] = str(k)+"\n"

        #/////creo il file txt per poterlo modificere come testo e poi passarlo in un json///
        open(namefile+'temp.json', 'w', encoding='utf-8')
        with open(namefile+'temp.json', 'a', encoding='utf-8') as f:
            # f.write("lol"+dizionariobello['nodes'][1]['data'][1]['label'])
            #f.write('{"elements":')
            #spaghetti=list(G.edges)
            f.write('{"elements":')
            f.write('{')
            f.write('"nodes": [')
            #indicenodotrovato
            for x in range(n_nodi):
                nodoattuale=lista_deinodi[x]
                f.write('{')
                f.write('"data":')
                f.write('  {')
                f.write('"id":"'+str(id_nodi_acasissimo[x])+'",')
                f.write('"label":"'+ lista_deinodi[x]+'",')
                f.write('"functann":"'+dizzi[nodoattuale][1]+'",')
                f.write('"nameEC":"'+dizzi[nodoattuale][2]+'",')
                f.write('"gn":"'+dizzi[nodoattuale][3]+'",')
                f.write('"vit3":"'+dizzi[nodoattuale][4]+'",')
                f.write('"nw1":"'+dizzi[nodoattuale][5]+'",')
                f.write('"nw2":"'+dizzi[nodoattuale][6]+'",')
                f.write('"nodec":"'+dict_isoColor[x+1]+'"')
                f.write('}')
                f.write('},')
            f.write('],')
            f.write('"edges": [')
            for x in range(n_Archi):
                f.write('{')
                f.write('"data":')
                f.write('  {')
                f.write('"id":"e'+ str(x)+'",')
                f.write('"source":"'+str(lista_degliarchi[x][0])+'",')
                f.write('"target":"'+str(lista_degliarchi[x][1])+'",')
                f.write('"weight":"'+'",')
                f.write('"codec":')

                if (spaghetti[x][0],spaghetti[x][1]) in estrongPos:
                    f.write('"1"')
                elif (spaghetti[x][0],spaghetti[x][1]) in estrongNeg:
                    f.write('"2"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumPos:
                    f.write('"3"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumNeg:
                    f.write('"4"')
                elif  (spaghetti[x][0],spaghetti[x][1])in eweakPos:
                    f.write('"5"')
                else:
                    f.write('"6"')

                f.write('}')
                f.write('},')
            f.write(']')
            f.write('}')
            f.write('}')
            f.close

            #da printare su terminale//////
            '''
            spaghetti=list(G.edges)
            print('{')
            print('\t'+'"nodes": [')
            for x in range(n_nodi):
                print('\t\t\t'+'   {')
                print('\t\t\t\t'+'"data":')
                print('\t\t\t\t'+'  {')
                print('\t\t\t\t\t'+'"id":"'+ str(x)+'",')
                print('\t\t\t\t\t'+'"label":"'+ lista_deinodi[x]+'",')
                print('\t\t\t\t\t'+'"vit3":"'+'",')
                print('\t\t\t\t\t'+'"vit2":"'+'"')
                print('\t\t\t\t'+'}')
                print('\t\t\t     },')
            print('\t\t\t],')
            print('\t'+'"edges": [')
            for x in range(n_Archi):
                print('\t\t\t'+'   {')
                print('\t\t\t\t'+'"data":')
                print('\t\t\t\t'+'  {')
                print('\t\t\t\t\t'+'"id":"'+ str(x)+'",')
                print('\t\t\t\t\t'+'"source":"'+lista_degliarchi[x][0]+'",')
                print('\t\t\t\t\t'+'"target":"'+lista_degliarchi[x][1]+'",')
                print('\t\t\t\t\t'+'"wheight":"'+'",')
                print('\t\t\t\t\t'+'"colourId":')
                if (spaghetti[x][0],spaghetti[x][1]) in estrongPos:
                    print('"estrongpos"')
                elif (spaghetti[x][0],spaghetti[x][1]) in estrongNeg:
                    print('"estrongneg"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumPos:
                    print('"emediumpos"')
                elif (spaghetti[x][0],spaghetti[x][1]) in emediumNeg:
                    print('"emediumneg"')
                elif  (spaghetti[x][0],spaghetti[x][1])in eweakPos:
                    print('"eweakpos"')
                else:
                    print('"eweakneg"')
                print('\t\t\t\t'+'}')
                print('\t\t\t     },')
            print('\t\t\t]')
            print('}')
            print(spaghetti)
            #print(lista_deinodi)
            #print(lista_degliarchi)
            print(estrongNeg)
            #print(ePos)
            #print(eNeg)
            '''

        #modifico json///////////////////////////////////
        polenta=[]
        with open(namefile+'temp.json', 'r', encoding='utf-8') as f:
            polenta=f.read()
            #print(polenta)
        polentaold=polenta
        pattern1= re.compile(r'\,\]')
        pattern2= re.compile(r'\,')
        pattern3=re.compile(r'\[')
        pattern4=re.compile(r'\]')
        pattern5=re.compile(r'\{')
        pattern6=re.compile(r'\}')
        pattern7=re.compile(r'\n\"')

        if pattern1.search(polenta):
            polenta=pattern1.sub(r']',polenta)
        if pattern7.search(polenta):
            polenta=pattern7.sub(r'"',polenta)
        if pattern2.search(polenta):
            polenta=pattern2.sub(r',\n',polenta)
        if pattern3.search(polenta):
            polenta=pattern3.sub(r'[\n',polenta)
        if pattern4.search(polenta):
            polenta=pattern4.sub(r']\n',polenta)
        if pattern5.search(polenta):
            polenta=pattern5.sub(r'{\n',polenta)
        if pattern6.search(polenta):
            polenta=pattern6.sub(r'}\n',polenta)


        stampare=polenta
        f=open(namefile+'final.json', 'w', encoding='utf-8')
        f.write(stampare)
        f.close
        #print(polenta)

        '''
        filetesto = open(namefile+"yo.txt","r")
        polentaold=filetesto.read()
        polenta=polentaold
        pattern= re.compile(r',]')
        with open("yo.txt") as f:
            for line in f:
                if pattern.search(line):
                    polenta=pattern.sub(r']',polenta)
        '''
        os.remove(namefile+'temp.json')
    #/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    #Clean graph and pyplot
    G.clear()
    plt.clf()
    plt.close()




#Draw Venn Diagram of common genes
def printVenn(listForVenn, couples, nameDir):
    for k in listForVenn:
        listKey = sorted([u for u in k.keys()], key=len)

        #Draw venn diagram with 2 sets
        if len(couples[listForVenn.index(k)]) == 2:
            elemSubSets = (len(k[str(listKey[0])]), len(k[str(listKey[1])]), int(len(k[str(listKey[2])])/2))
            v = venn2(subsets = elemSubSets, set_labels = ((re.findall(r'\w+', str(listKey[0])))[0], (re.findall(r'\w+', str(listKey[1])))[0]))
            plt.savefig(nameDir+str(listForVenn.index(k))+'/venn.png')
            print('Create: \''+str(listForVenn.index(k))+'/venn.png\'')
        #Draw venn diagram with 3 sets
        if len(couples[listForVenn.index(k)]) == 3:
            elemSubSets = ('',)
            if len(listKey) == 7:
                listKey = sorted(listKey[:3])+sorted(listKey[3:-1], key=itemgetter(0,1))+[listKey[-1]]
                #(A,B,AB,C,AC,BC,ABC)
                elemSubSets = (len(k[str(listKey[0])]), len(k[str(listKey[1])]), int(len(k[str(listKey[3])])/2), len(k[str(listKey[2])]), int(len(k[str(listKey[4])])/2), int(len(k[str(listKey[5])])/2), int(len(k[str(listKey[6])])/3))
                v = venn3(subsets = elemSubSets, set_labels = ((re.findall(r'\w+', str(listKey[0])))[0], (re.findall(r'\w+', str(listKey[1])))[0], (re.findall(r'\w+', str(listKey[2])))[0]))
            else:
                pass
                #TODO: find solution if there isn't some element in listKey
            plt.savefig(nameDir+str(listForVenn.index(k))+'/venn.png')
            print('Create: \''+str(listForVenn.index(k))+'/venn.png\'')
        #Draw Venn diagram with 4 sets
        if len(couples[listForVenn.index(k)]) == 4:
            #made all possible combination of genes
            listKey = sorted([str([elem]) for elem in couples[listForVenn.index(k)]])
            listKeyL = sorted([[elem] for elem in couples[listForVenn.index(k)]])
            listKey.append(str([listKeyL[0][0],listKeyL[1][0],listKeyL[2][0],listKeyL[3][0]]))
            i = 0
            j = 1
            z = 2
            while i < 4:
                tmp = listKeyL[i][0]
                while j < 4:
                    listKey.append(str([tmp,listKeyL[j][0]]))
                    while z < 4:
                        listKey.append(str([tmp,listKeyL[j][0],listKeyL[z][0]]))
                        z+=1
                    j+=1
                    z=j+1
                i+=1
                j=i+1
                z=j+1
            listKey = sorted(listKey, key=len)
            keyAndNumber = [(key,len(k[key])) for key in sorted(k.keys(), key=len)]
            completeKeyNumber = {}
            for key in listKey:
                if key in [elem[0] for elem in keyAndNumber]:
                    completeKeyNumber[key] = keyAndNumber[[elem[0] for elem in keyAndNumber].index(key)][1]
                else:
                    completeKeyNumber[key] = 0
            #print(completeKeyNumber)
            dictLabels = {
                '0001': completeKeyNumber[listKey[0]],
                '0010': completeKeyNumber[listKey[1]],
                '0011': int(completeKeyNumber[listKey[4]]/2),
                '0100': completeKeyNumber[listKey[2]],
                '0101': int(completeKeyNumber[listKey[5]]/2),
                '0110': int(completeKeyNumber[listKey[7]]/2),
                '0111': int(completeKeyNumber[listKey[10]]/3),
                '1000': completeKeyNumber[listKey[3]],
                '1001': int(completeKeyNumber[listKey[6]]/2),
                '1010': int(completeKeyNumber[listKey[8]]/2),
                '1011': int(completeKeyNumber[listKey[11]]/3),
                '1100': int(completeKeyNumber[listKey[9]]/2),
                '1101': int(completeKeyNumber[listKey[12]]/3),
                '1110': int(completeKeyNumber[listKey[13]]/3),
                '1111': int(completeKeyNumber[listKey[14]]/4),
            }
            v = vennD.venn4(dictLabels, names=[listKeyL[3][0],listKeyL[2][0],listKeyL[1][0],listKeyL[0][0]])
            plt.savefig(nameDir+str(listForVenn.index(k))+'/venn.png')
            print('Create: \''+str(listForVenn.index(k))+'/venn.png\'')

        #Draw Venn diagram with 5 sets
        if len(couples[listForVenn.index(k)]) == 5:
            #made all possible combination of genes
            listKey = sorted([str([elem]) for elem in couples[listForVenn.index(k)]])
            listKeyL = sorted([[elem] for elem in couples[listForVenn.index(k)]])
            listKey.append(str([listKeyL[0][0],listKeyL[1][0],listKeyL[2][0],listKeyL[3][0],listKeyL[4][0]]))
            i = 0
            j = 1
            w = 3
            z = 2
            while i < 5:
                tmp = listKeyL[i][0]
                while j < 5:
                    listKey.append(str([tmp,listKeyL[j][0]]))
                    while z < 5:
                        listKey.append(str([tmp,listKeyL[j][0],listKeyL[z][0]]))
                        while w < 5:
                            listKey.append(str([tmp,listKeyL[j][0],listKeyL[z][0],listKeyL[w][0]]))
                            w+=1
                        z+=1
                        w=z+1
                    j+=1
                    z=j+1
                    w=z+1
                i+=1
                j=i+1
                z=j+1
                w=z+1
            listKey = sorted(listKey, key=len)
            keyAndNumber = [(key,len(k[key])) for key in sorted(k.keys(), key=len)]
            completeKeyNumber = {}
            for key in listKey:
                if key in [elem[0] for elem in keyAndNumber]:
                    completeKeyNumber[key] = keyAndNumber[[elem[0] for elem in keyAndNumber].index(key)][1]
                else:
                    completeKeyNumber[key] = 0
            #print(completeKeyNumber)
            dictLabels = {
                '00001': completeKeyNumber[listKey[0]],
                '00010': completeKeyNumber[listKey[1]],
                '00011': int(completeKeyNumber[listKey[5]]/2),
                '00100': completeKeyNumber[listKey[2]],
                '00101': int(completeKeyNumber[listKey[6]]/2),
                '00110': int(completeKeyNumber[listKey[9]]/2),
                '00111': int(completeKeyNumber[listKey[15]]/3),
                '01000': completeKeyNumber[listKey[3]],
                '01001': int(completeKeyNumber[listKey[7]]/2),
                '01010': int(completeKeyNumber[listKey[10]]/2),
                '01011': int(completeKeyNumber[listKey[16]]/3),
                '01100': int(completeKeyNumber[listKey[12]]/2),
                '01101': int(completeKeyNumber[listKey[18]]/3),
                '01110': int(completeKeyNumber[listKey[21]]/3),
                '01111': int(completeKeyNumber[listKey[25]]/4),
                '10000': completeKeyNumber[listKey[4]],
                '10001': int(completeKeyNumber[listKey[8]]/2),
                '10010': int(completeKeyNumber[listKey[11]]/2),
                '10011': int(completeKeyNumber[listKey[17]]/3),
                '10100': int(completeKeyNumber[listKey[13]]/2),
                '10101': int(completeKeyNumber[listKey[19]]/3),
                '10110': int(completeKeyNumber[listKey[22]]/3),
                '10111': int(completeKeyNumber[listKey[26]]/4),
                '11000': int(completeKeyNumber[listKey[14]]/2),
                '11001': int(completeKeyNumber[listKey[20]]/3),
                '11010': int(completeKeyNumber[listKey[23]]/3),
                '11011': int(completeKeyNumber[listKey[27]]/4),
                '11100': int(completeKeyNumber[listKey[24]]/3),
                '11101': int(completeKeyNumber[listKey[28]]/4),
                '11110': int(completeKeyNumber[listKey[29]]/4),
                '11111': int(completeKeyNumber[listKey[30]]/5),
            }
            v = vennD.venn5(dictLabels, names=[listKeyL[4][0],listKeyL[3][0],listKeyL[2][0],listKeyL[1][0],listKeyL[0][0]])
            plt.savefig(nameDir+str(listForVenn.index(k))+'/venn.png')
            print('Create: \''+str(listForVenn.index(k))+'/venn.png\'')
        plt.clf()
        plt.close()

#Draw histogram graph of common genes
def printHistogram(listCommonGenes, listFiles, nameDir):
    for l in listCommonGenes:
        listFrel = []
        for k in l[0]:
            listFrel.append([k])
        i = 1
        while i < len(l):
            listFrel[[u[0] for u in listFrel].index(l[i][0])].append(l[i][3])
            i += 1
        frelOriginalFile = []
        for name in [u[0] for u in listFiles]:
            if name in [u[0] for u in listFrel]:
                tmp = [name]
                for u in listFiles:
                    if u[0] == name:
                        tmp = tmp + [k[2] for k in u[1:]]
                frelOriginalFile.append(tmp)
        listFrel = sorted(listFrel, key=ut.ord)
        frelOriginalFile = sorted(frelOriginalFile, key=ut.ord)
        #Prepare parameters to draw the histograms
        num_rows = int((float(len(listFrel))-0.1)/2)+1
        if num_rows == 1:
            num_rows+=1;
        num_bins = 20
        fig, axes = plt.subplots(num_rows, 2, tight_layout=True)

        #Draw each histograms
        counter = 0
        for i in range(num_rows):
            for j in range(2):
                ax = axes[i][j]
                # Plot when we have data
                if counter < len(listFrel):
                    #Draw histogram with linear axes
                    ax.hist([frelOriginalFile[counter][1:], listFrel[counter][1:]], bins=num_bins, range=[0.0, 1.0], edgecolor='black', linewidth=1.2, color=['green', 'blue'], alpha=0.5, label=[str(listFrel[counter][0])+' original list', str(listFrel[counter][0])+' analyzed genes'])
                    # ax.hist(frelOriginalFile[counter][1:], bins=num_bins, density=True, edgecolor='black', linewidth=1.2, color='green', alpha=0.5, label='{}'.format(listFrel[counter][0]))
                    # ax.hist(listFrel[counter][1:], bins=num_bins, density=True, edgecolor='black', linewidth=1.2, color='blue', alpha=0.5, label='{}'.format(listFrel[counter][0]))
                    ax.set_xlabel('Frequency')
                    ax.set_ylabel('Number genes')
                    leg = ax.legend(loc='upper left')
                    leg.draw_frame(False)
                # Remove axis when we no longer have data
                else:
                    ax.set_axis_off()
                counter += 1

        #create dir for each couple of genes
        nameDirGenes = nameDir+str(listCommonGenes.index(l))+'/'
        # tmp = plt.gcf().get_size_inches()
        plt.gcf().set_size_inches(15, 10)
        plt.savefig(nameDirGenes+'histogram.png')
        print('Create: \''+nameDirGenes+'histogram.png\'')
        # plt.gcf().set_size_inches(tmp)
        plt.clf()
        plt.close()


def printDegree(G, idNode, namefile):
    f = open(namefile.split("Graph")[0]+'degree_nodes.txt', 'w')
    f.write('Node,Degree\n')
    for key in idNode.keys():
        f.write(key+','+str(G.degree[idNode[key]])+'\n')
    f.close()


'''
def nodefounder(a,b,c):
    quantograndelistlinename=len(a)
    for y in range(quantograndelistlinename):
        nododatrovare=b[c]
        indicetrovato=y
        #global indicenodotrovato
        if (nododatrovare==a[y][0]):
            #indicetrovato=y
            print(indicetrovato)
            return indicetrovato'''
