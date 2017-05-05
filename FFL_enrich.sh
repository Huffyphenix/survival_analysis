#!/usr/bin/sh
DATA_PATH=`pwd`
id=$1
#for id in `cat /project/huff/huff/TCGA_survival/data/cancer.list.txt`
#do
mkdir $DATA_PATH/$id/FFL_filter
mkdir $DATA_PATH/$id/enrich
cat $DATA_PATH/$id/${id}_miRNA_survival.p_signif.info.txt|cut -f 1 |cut -d '|' -f 1 >$DATA_PATH/$id/FFL_filter/mirna.txt
cat $DATA_PATH/$id/${id}_RNA_survival.p_signif.info.txt |cut -f 1 |cut -d '|' -f 1 >$DATA_PATH/$id/${id}_RNA_survival.p_signif.list
awk -F'\t' 'NR==FNR{a[$1]=$1}NR>FNR{print a[$1]}' /project/huff/huff/id_corresponding/AnimalTFDB_id/AnimalTFDB_tf $DATA_PATH/$id/${id}_RNA_survival.p_signif.list|grep -v '^$' >$DATA_PATH/$id/FFL_filter/TF.txt
grep -Fw -v -f $DATA_PATH/$id/FFL_filter/TF.txt $DATA_PATH/$id/${id}_RNA_survival.p_signif.list >$DATA_PATH/$id/FFL_filter/Protein.txt
bash /project/huff/huff/script/scripts_for_network/construct_FFL_network.sh -m $DATA_PATH/$id/FFL_filter/mirna.txt -t $DATA_PATH/$id/FFL_filter/TF.txt -g $DATA_PATH/$id/FFL_filter/Protein.txt -a /project/huff/huff/data/FFL_database/TF2PrGene/TF2gene.all -b /project/huff/huff/data/FFL_database/TF2miRNA/TF2miRNA.all -c /project/huff/huff/data/FFL_database/miRNA2gene/miRNA2gene.all -d $DATA_PATH/$id/FFL_filter
cat $DATA_PATH/$id/${id}_RNA_survival.p_signif_exp1000_FC2.info.txt |cut -f 1|cut -d '|' -f 1 >$DATA_PATH/$id/${id}_RNA_survival.p_signif.exp1000_FC2.list.symbol
python /project/zhangq/tissue_special/analysis/scripts/hyper_test.py -e $DATA_PATH/$id/${id}_RNA_survival.p_signif.exp1000_FC2.list.symbol -b /project/zhangq/tissue_special/analysis/scripts/test/human/KO_human_zq_20161216kegg_version.desc -c 3 -f 0.1 -o $DATA_PATH/$id/enrich/${id}_RNA_survival.p_signif.exp1000_FC2.enrich
cat $DATA_PATH/$id/${id}_RNA_survival.p_signif.list.txt |cut -d '|' -f 1  >$DATA_PATH/$id/${id}_RNA_survival.p_signif.list.symbol
python /project/zhangq/tissue_special/analysis/scripts/hyper_test.py -e $DATA_PATH/$id/${id}_RNA_survival.p_signif.list.symbol -b /project/zhangq/tissue_special/analysis/scripts/test/human/KO_human_zq_20161216kegg_version.desc -c 3 -f 0.1 -o $DATA_PATH/$id/enrich/${id}_RNA_survival.p_signif.enrich
Rscript /project/huff/huff/TCGA_survival/script/target.R $DATA_PATH/$id/FFL_filter/miRNA2gene $DATA_PATH/$id/FFL_filter/TF2gene $DATA_PATH/$id/FFL_filter/miRNA2TF $DATA_PATH/$id/FFL_filter/TF2miRNA $DATA_PATH/$id/${id}_RNA_survival.p_signif.info.txt $DATA_PATH $id $DATA_PATH/$id/${id}_miRNA_survival.p_signif.info.txt
#done


