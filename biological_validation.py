import sys
import os
import re
import datetime as dt
import rpy2.robjects as ro
import pandas as pd
import lib.diffexp_go_analysis as topGO

#Print info about cmd command if the call is wrong
def printInfo():
    print('Usage: python3 biological_validation.py PARAM LIST_GENES')
    print('\nPARAM:')
    print('\t-topGO\t\tExecute GO validation')
    print('\t-xstreme\t\tPrepare file for XSTREME analysis')
    print('\nLIST_GENES:\tList of genes. Esample: An example is the file: \'example_lists/default_example_STS/STS_biologicalValidation.csv\'')
    print('\n-->TEST EXAMPLE TOPGO: python3 biological_validation.py -topGO example_lists/default_example_STS/STS_biologicalValidation.csv')
    print('\n-->TEST EXAMPLE XSTREME: python3 biological_validation.py -xstreme example_lists/default_example_STS/STS_biologicalValidation.csv')
    sys.exit(-1)

#Create saving directory
def createSavingDir():
    #create directory
    if not os.path.exists('outputBiologicalValidation'):
        os.mkdir('outputBiologicalValidation')
    #write graph in a .txt file
    time = dt.datetime.now()
    printTime = time.strftime("%Y%m%d%H%M%S")
    os.mkdir('outputBiologicalValidation/'+printTime)
    global nameDir
    nameDir = 'outputBiologicalValidation/'+printTime+'/'
    print('Creating directory: \''+nameDir+'\'', flush=True)
    return nameDir

#Print txt with all information of topGO
def print_output_topGO(results_table, results, nameDir, ontology):
    #Write file .txt
    f = open(nameDir+'validation_topGO_'+ontology+'.txt', 'w')
    f.write(str(results)+'\n\n')
    pd_dt = ro.conversion.rpy2py(results_table)
    f.write(str(pd_dt))
    f.close();
    print('Create: '+nameDir+'validation_topGO_'+ontology+'.txt', flush=True)

#Create fasta file
def createFasta(listGenes, completeFasta, dir):
    with open(listGenes) as in_handle:
        genes_list = topGO.parse_input_csv(in_handle)
    genes_list=[y for k in [re.split('<BR>',g) for g in genes_list] for y in k]
    f_complete = open(completeFasta, 'r')
    completeF = f_complete.read().split('\n')
    dict_FASTA = {}
    i = 0
    while i < len(completeF):
        if completeF[i][1:].upper() in genes_list:
            dict_FASTA[completeF[i]] = completeF[i+1]
        i += 2
    f_complete.close()
    namefileFASTA = nameDir+'fasta_'+listGenes.split('/')[-1][:-4]+'.fasta'
    print(namefileFASTA)
    f = open(namefileFASTA, 'w')
    for k in dict_FASTA.keys():
        f.write(k+'\n'+dict_FASTA[k]+'\n')
    f.close()
    return namefileFASTA

def main():
    if len(sys.argv) >= 1:
        if sys.argv[1] == '-topGO':
            #create saving directory
            nameDir = createSavingDir()
            #results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt')
            #Run BP ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt', nameDir, 'BP')
            print_output_topGO(results_table, results, nameDir, 'BP')
            #Run MF ontology validation
            results_table, results = topGO.topGO_analysis(sys.argv[2], 'import_doc/V1_GOcomplete.txt', nameDir, 'MF')
            print_output_topGO(results_table, results, nameDir, 'MF')
        elif sys.argv[1] == '-xstreme':
            #create saving directory
            nameDir = createSavingDir()
            #create fasta file from list of genes
            fastaFile = createFasta(sys.argv[2], 'import_doc/grape_1k_upstream.fasta', nameDir)
        elif sys.argv[1] == '--help':
            printInfo()
        else:
            print('ERROR PARAM')
            printInfo()
    else:
        print('ERROR: wrong nuomber of parameters')

#Calls the main() function.
if __name__ == '__main__':
    main()
