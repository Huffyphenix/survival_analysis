for id in `cat /project/huff/huff/TCGA_survival/data/cancer.list.txt`
do
nohup bash /project/huff/huff/TCGA_survival/script/survival.sh -e /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_$id.mRNAseq_Preprocess.Level_3.2016012800.0.0/$id.uncv2.mRNAseq_RSEM_all.txt -s /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_$id.Merge_Clinical.Level_1.2016012800.0.0/$id.clin.merged.txt.2 -c $id -l RNA -p 0.05 -b 0.5  > $id.RNA.out 2>&1 &
done

