# This script downloads compressed 'tree sequence' representations
# of the public 1000 genomes sequence dataset. 
# These files will take up 2GB of hard disk space
#
# More information about these tree sequences can be found at:
# https://zenodo.org/record/3051855#.YMzxeZOpE5Y
# https://tskit.dev/tskit/docs/stable/python-api.html

wget https://zenodo.org/record/3051855/files/1kg_chr1.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr2.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr3.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr4.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr5.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr6.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr7.trees.tsz 
wget https://zenodo.org/record/3051855/files/1kg_chr8.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr9.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr10.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr11.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr12.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr13.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr14.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr15.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr16.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr17.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr18.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr19.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr20.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr21.trees.tsz
wget https://zenodo.org/record/3051855/files/1kg_chr22.trees.tsz
mv *.tsz data
