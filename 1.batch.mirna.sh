for id in `cat /project/huff/huff/TCGA_survival/data/cancer.list.txt`
do
nohup bash /project/huff/huff/TCGA_survival/script/survival.sh -e /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_$id.miRseq_Mature_Preprocess.Level_3.2016012800.0.0/$id.miRseq_mature_RPM.txt -s /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_$id.Merge_Clinical.Level_1.2016012800.0.0/$id.clin.merged.txt.2 -c $id -l miRNA -p 0.05 -b 0.5  > $id.miRNA.out 2>&1 &
done

id=LAML
nohup bash /project/huff/huff/TCGA_survival/script/survival.sh -e /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_LAML.miRseq_Preprocess.Level_3.2016012800.0.0/LAML.miRseq_RPKM.txt -s /project/huff/huff/TCGA_survival/data/$id/gdac.broadinstitute.org_$id.Merge_Clinical.Level_1.2016012800.0.0/$id.clin.merged.txt.2 -c $id -l miRNA -p 0.05 -b 0.5  > $id.miRNA.out 2>&1 &
