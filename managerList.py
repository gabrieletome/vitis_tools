import sys
import lib.utilities as ut
import lib.vitis as vitis

def main():
    #Get the parameters from the command line, using printInfo() as a fallback
    #Minimun 2 parameters: [0] = name file python, [1] = name file of the list of the gene (.csv/.zip)
    if len(sys.argv) >= 2:
        if sys.argv[1] == '--help':
            ut.printInfo()
        elif sys.argv[1] == '-vitis':
            #read parameters
            cmd = ut.readParameters(sys.argv)
            #build matrix of genes
            matrixGenes = vitis.buildMatrixGenesVitis(cmd[0], cmd[1])
            print('Build matrix successfully', flush=True)
            #build graph of interaction between genes
            print("Building complete graph...", flush=True)
            graphGenes = ut.manageDuplicates(ut.buildGraph(matrixGenes))
            print('Build graph successufully', flush=True)
            #find core of graphGenes
            print("Searching core network...", flush=True)
            coreGraph = ut.findCoreGraph(graphGenes)
            print('Core find successufully', flush=True)
            #print in output complete graph and core
            vitis.printOutput(coreGraph, graphGenes)
        else:
            print('ERROR: Need 1 parameters')
            print('Usage: python3 managerList.py PARAM [FILTERS]... -files [FILES]...')
            print('PARAM:')
            print('\t-vitis\tLists of vitis')
    else:
        print('ERROR: wrong format')
        ut.printInfo()

#Calls the main() function.
if __name__ == '__main__':
    main()
