#install python libraries
pip3 install -r import_doc/requirements.txt
#install BiocManager and topGO
sudo R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); library('BiocManager'); BiocManager::install('topGO')"
