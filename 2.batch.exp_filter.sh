for id in `cat /project/huff/huff/TCGA_survival/data/cancer.list.txt`
do
nohup bash /project/huff/huff/TCGA_survival/script/FFL_enrich.sh $id >$id.FFL_enrich.out 2>&1 &
done

