# Tools_Expand_Genes_Network
 Tools used to analyze and expand gene networks of Vitis Vinifera. All this tools used Python 3.

## Prerequisites
 Python 3
 Use requirements.txt to install the required library:
 ```
  pip3 install -r requirements.txt
 ```
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


 For *biological_validation.py* is required R with [*topGO*](https://bioconductor.org/packages/release/bioc/html/topGO.html) and [DREME](http://meme-suite.org/doc/dreme.html) form MEME suite.
   Install *topGO* library:
   ```
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

    BiocManager::install("topGO")
   ```
   Install DREME: [DREME tool](http://meme-suite.org/doc/download.html)

## TOOL 1: *managerList.py*
Tool used to build the complete graph of interaction and find the LGN.
How to use:
```
Usage: python3 managerList.py -vitis [FILTERS]... -files [FILES]...
FILES can be a list of .csv or .zip. These are the expansion list from OneGenE
FILTERS:
 '-a'             Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS
 '-f [NUMBER]'    Ignored genes with frel<=NUMBER
 '-t [PATTERN,...]'      Take genes that in 'Network1' or 'Network2' column there is at least one pattern
```
Example with Cuticle genes:
```
python3 managerList.py -vitis -f 0.1 -a -files example_lists/Vitis_7genes_MYB-ERF/Example.zip
```
Help:
```
python3 managerList.py --help
```

## TOOL 2: *integrateCoupleGenes.py*
Tool used to expand the LGN.
How to use:
```
Usage: python3 integrateCoupleGenes.py -vitis TYPEA [FILTERS]... -files [GENES] [FILES]...
TYPEA:
 '-frel'                Build expansion network based on FREL. Required filter '-f'
 '-rank [INT]'          Build expansion network based on RANK. Take top genes
 '-shared'              Build expansion network based on SHARED GENES
 '-pattern [PATTERNS]'  Build expansion network based on genes that have at least one patterns in 'Network1' or 'Network2' columns
FILTERS:
 '-a'           Autosave image of graphs. If -a is present, it save automatically .png. USE IN MICROSOFT WINDOWS
 '-c'           Add edges between associated genes
 '-e'           Print Venn Diagram and Histogram for complete analysis
 '-f [NUMBER]'  Ignored genes with frel<=NUMBER
GENES: file .csv with the genes to analyze. Example: ('CoupleGeneToIntegrate/coupleGene.csv')
FILES can be a list of .csv or .zip
```
Example with Cuticle genes:
```
python3 integrateCoupleGenes.py -vitis -shared -f 0.1 -e -a -files CoupleGeneToIntegrate/coupleGene0.csv example_lists/Vitis_7genes_MYB-ERF/Example.zip

python3 integrateCoupleGenes.py -vitis -pattern vv60 -f 0.1 -e -a -files CoupleGeneToIntegrate/GroupB_STS.csv example_lists/Esempio_STS/StSyMYB.zip
```
Help:
```
python3 integrateCoupleGenes.py --help
```

## TOOL 3: *biological_validation.py*
Tool that execute biological validation.
How to use:
```
Usage: python3 biological_validation.py PARAM [FILTERS]... LIST_GENES COMPLETE_GENOME
PARAM:
 '-topGO' Execute GO validation
 '-dreme' Execute DREME analysis
FILTERS (ONLY FOR DREME):
 '-exe'  Execute dreme analysis in local. Without prepare fasta file for DREME website
LIST_GENES      List of genes in .csv file.
COMPLETE_GENOME COmplete file information for validation
 * For topGO: file map from vitis ID to GO ID. ('import_doc/V1_GOcomplete.txt')
 * For Dreme: complete list of genes in genome. ('import_doc/grape_1k_upstream.fasta')
```
Example topGO with Cuticle genes:
```
python3 biological_validation.py -topGO example_lists/Vitis_7genes_MYB-ERF/topGO_Vitis/MYB_ERF_topGO0.csv import_doc/V1_GOcomplete.txt
```
Example DREME with Cuticle genes:
```
python3 biological_validation.py -dreme example_lists/Vitis_7genes_MYB-ERF/topGO_Vitis/MYB_ERF_topGO0.csv import_doc/grape_1k_upstream.fasta
```
Help:
```
python3 biological_validation.py --help
```
