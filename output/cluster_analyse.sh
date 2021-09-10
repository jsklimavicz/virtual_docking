#! /bin/bash

echo "$(date) Commencing cluster analysis." >> /tmp/analysis.log
rm res_clust*.pdb lig_clust*.mol2 ligands_cluster*.csv

echo Preparing pca parameters and performing Bayesian clustering...
python3 -u ./bayesMM2.py output_filtered.csv -i 10 -v -V 0.90 -c 24 -m ward
echo Done!

echo Graphing clusters...
Rscript graph_clusters.R 
echo Done!

echo Preparing cluster pdb files...
for file in ligands_cluster*.csv; do
#    echo "Input: $file"
    cluster=$(ls $file | sed -e s/[^0-9]//g)
    python3 batch_pdb_merge.py $file DmelTyrR2_rigid.pdbqt -o res_clust$cluster lig_clust$cluster -n 25
done

#tar -czf clusters.tar.gz res_clust*.pdb lig_clust*.mol2 ligands_cluster*.csv DmelTyrR2_rigid.pdbqt cluster_graph*.png
#echo Done! Files zipped as clusters.tar.gz

echo "$(date) Cluster analysis successful." >> /tmp/analysis.log

echo "Making grouped compound PDF. This may take a while..."
python3 -u latex_writer.py bayesMM_ligands.csv -p 10 -v
echo PDF Completed.
echo "$(date) PDF write successful." >> /tmp/analysis.log

tar -czf clusters.tar.gz res_clust*.pdb lig_clust*.mol2 ligands_cluster*.csv DmelTyrR2_rigid.pdbqt cluster_graph*.png compound_structures.pdf
echo Done! Files zipped as clusters.tar.gz
