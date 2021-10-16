# Tools_Expand_Genes_Network
 Tools used to analyze and expand gene networks of Vitis Vinifera. All this tools use Python 3 and R.

## Prerequisites
 Python 3

 To install the required Python libraries and R packages:
 ```
  ./setup.sh
 ```

 In directory 'import_doc' is present the file 'vv_exprdata_2.csv'.  It is a large file uploaded with Git LFS that generate a pointer to the real file. To have the correct version of the file you need to:
  - download the repository with Git LFS, or,
  - get the file at this link (<!--insert link-->) and substitute to the file in 'import_doc' obtained with git clone.
 <!-- Use *requirements.txt* to install the required libraries:
 ```
  pip3 install -r import_doc/requirements.txt
 ``` -->
 <!-- Library of Python (install with pip3):
   * *datetime*
   ```
    pip3 install datetime
   ```
   * *pandas*
   ```
    pip3 install pandas
   ```
   * *rpy2*
   ```
    pip3 install rpy2
   ```
   * *matplotlib*
   ```
    pip3 install matplotlib
   ```
   * *matplotlib-venn*
   ```
    pip3 install matplotlib-venn
   ```
   * *networkx*
   ```
    pip3 install networkx
   ```
   * *numpy*
   ```
    pip3 install numpy
   ```
   * *scipy*
   ```
    pip3 install scipy
   ```
   * *rpack*
   ```
    pip3 install rectangle-packer
   ``` -->

   <!-- Install *topGO* library:
   ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("topGO")
   ``` -->

## TOOL 1: *managerList.py*
Tool used to build the complete graph of interaction and find the LGN.
The output of this tool is the graphical representation (*.png*) of the complete graph discovered observing the file in input with corresponding legend (*.txt*). It is created the *.json* file compatible with Cytoscape. The tool generates a *.csv* file containing the textual representation of the graph, basically are written the edges of the graph. ManagerList returns the degree of the genes and the connected components of the graph.
The LGN, Local Gene Network of the previous graph, discovered using an adapted version of Charikar algorithm, corresponds to the *Core graph* files. The files for the LGN are the same of the complete graph previously described.
How to use:
```
Usage: python3 managerList.py -vitis [FILTERS]... -files [FILES]...

FILTERS:
 '-a'             Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS WSL
 '-f [NUMBER]'    Ignored genes with frel<=NUMBER
 '-t [PATTERN,...]'      Take genes that in 'Network1' or 'Network2' column there is at least one pattern

FILES: can be a list of .csv or .zip. These are the expansion list from OneGenE

```
Example with STS:
```
python3 managerList.py -vitis -f 0.5 -a -files example_lists/default_example_STS/STS_example.zip
```
Help:
```
python3 managerList.py --help
```

## TOOL 2: *integrateCoupleGenes.py*
Tool used to expand the LGN adding to the graph the genes connected to the LGN with specific characteristic. The expansion can be done by:
  - apply a filter on the relative frequency;
  - apply a filter on the rank of the genes in the expansion list;
  - taking only the common genes;
  - apply a filter on the columns *Network1* and *Network2*
In output, the tool returns the graphical (*graphGenes.png*) and textual (*edges_graph.csv*) representation of the expanded graph  with the corresponding legend and the file *.json* to import in Cytoscape. Optionally, it can be obtained the histogram and the venn diagram. The histogram represents the genes selected for the expansion (blue columns) respect all the genes present in the expansion list of the starting gene (green columns). The venn diagram (graphical representation up to 5 genes in LGN, textual for any number of genes in LGN) represents how the genes in expansion list are shared between the LGN.
How to use:
```
Usage: python3 integrateCoupleGenes.py -vitis TYPEA [FILTERS]... -files [GENES] [FILES]...

TYPEA:
 '-frel'                Build expansion network applying a filter to the relative frequency of the edges. Required filter '-f'
 '-rank [INT]'          Build expansion network based on RANK of the genes in each expansion list given in input.
 '-shared'              Build expansion network based on SHARED GENES of the genes of the LGN.
 '-pattern [PATTERNS]'  Build expansion network based on genes that have at least one patterns in 'Network1' or 'Network2' columns

FILTERS:
 '-a'           Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS WSL
 '-c'           Add edges between associated genes
 '-e'           Print Venn Diagram and Histogram for complete analysis
 '-f [NUMBER]'  Ignored genes with frel<=NUMBER

GENES: file .csv with the genes to analyze. Example: ('CoupleGeneToIntegrate/default_example_STS.csv')

FILES: can be a list of .csv or .zip. These are the expansion list from OneGenE
```
Example with STS:
```
python3 integrateCoupleGenes.py -vitis -shared -f 0.1 -a -e -files CoupleGeneToIntegrate/default_example_STS.csv example_lists/default_example_STS/STS_example.zip
```
Help:
```
python3 integrateCoupleGenes.py --help
```

## TOOL 3: *biological_validation.py*
Tool that execute biological validation. It perform two analysis:
  - Computes a Gene Ontology analysis on set of genes.
  - Prepare the fasta file for the execution of XSTREME tool.
How to use:
```
Usage: python3 biological_validation.py PARAM LIST_GENES

PARAM:
 '-topGO' Execute GO validation
 '-xstreme' Execute XSTREME analysis

LIST_GENES      List of genes in .csv file. An example is the file: 'example_lists/default_example_STS/STS_biologicalValidation.csv'
```
Example Gene Ontology analysis:
```
python3 biological_validation.py -topGO example_lists/default_example_STS/STS_biologicalValidation.csv
```
Example command need for the preparation of files for XSTREME:
```
python3 biological_validation.py -xstreme example_lists/default_example_STS/STS_biologicalValidation.csv
```
Help:
```
python3 biological_validation.py --help
```
