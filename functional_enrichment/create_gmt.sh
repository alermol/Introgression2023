#!/bin/bash

# Script creates three GMT (Gene Matrix Transposed) files. For each GO category
# seperate file will be created.

mkdir temp/

echo "Downloading necessary files..."
wget -nv http://purl.obolibrary.org/obo/go.obo -P temp/
wget -nv https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/iwgsc_refseqv1.0_FunctionalAnnotation_v1.zip -P temp/

echo "Unzipping annotaion..."
cd temp/
unzip -qq iwgsc_refseqv1.0_FunctionalAnnotation_v1.zip
mv iwgsc_refseqv1.0_FunctionalAnnotation_v1/iwgsc_refseqv1.0_FunctionalAnnotation_v1__HCgenes_v1.0-repr.TEcleaned.TAB annotation.tab

echo "Creating GMT files..."
paste <(grep -o '^id: GO:[0-9]*' go.obo | sed 's/id:\ //') \
    <(grep -o '^name: .*' go.obo | sed 's/name:\ //') \
    <(grep -o '^namespace: [a-z_]*' go.obo | sed 's/namespace:\ //') > go_extracted.tsv

grep -o 'GO:[0-9]*' annotation.tab | sort -u > annot_go_terms.tab

while read line
do
    descr=$(grep $line go_extracted.tsv | cut -f 2 | tr ' ' '_')
    grep $line annotation.tab | cut -f 1 | sed 's/\.[1-9]$//g' | tr '\n' '\t' | awk -v l=$line -v d=$descr -v OFS='\t' '{print l,d,$0}'
done < annot_go_terms.tab  | tr '_' ' ' > result.gmt

grep 'biological_process' go_extracted.tsv | cut -f 1 | grep -f - result.gmt > ../GO:BP.gmt
grep 'molecular_function' go_extracted.tsv | cut -f 1 | grep -f - result.gmt > ../GO:MF.gmt
grep 'cellular_component' go_extracted.tsv | cut -f 1 | grep -f - result.gmt > ../GO:CC.gmt

mv result.gmt ../all_GO_categories.gmt

echo "GMT files were created..."
echo "Cleaning up..."

cd ..

rm -rf temp/

echo "Done"
