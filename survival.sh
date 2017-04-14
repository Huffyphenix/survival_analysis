#!/bin/sh
OPTIND=1
SCRIPTS_HOME=/project/huff/huff/TCGA_survival/script
BASE_PATH=$PWD
MLM_PVALUE=0.05
SUR_PVALUE=0.05
BOUND=0.5
usage() {
    cat <<EOF
Usage:bash $0 -e <exp_matrix> -s <TCGA_clinical_file> -c <cancer tyype> -l <label> -p [survival pvalue cutoff] -b [gene expression boundary] -d [workdir]
        -e <string> Expression matrix of the genes, colnames are gene id, rownames are sample id. Seperated by tab.
                    See input_example/survival.exp.
        -s <string> Raw TCGA clinical file: include survival time and status of the samples.
                    You don't need to handle it, just ensure it is in a proper format, see input_example/survival_time.txt.
        -c <cancer type> cancer type/abbreviation of cancer.
        -l <string> miR or RNA, or some label you want to show at the beginning your result filename.
        -p [numeric] Cutoff pvalue of survival analysis. Default is 0.05.
        -b [numeric] >0 and <1. The boundary of one gene's expression which you want to do survival analysis with. Default: 0.5.
        -d [string] working direcotry. Default is $PWD.
EOF
    exit 1
}
while getopts :e:s:c:l:p:b:d: ARGS; do
    case $ARGS in
                e) EXP_PATH=$OPTARG ;;
                s) SUR_TIME_PATH=$OPTARG ;;
                c) CANCER=$OPTARG ;;
                l) LABEL=$OPTARG ;;
                p) SUR_PVALUE=$OPTARG ;;
                b) BOUND=$OPTARG ;;
                d) BASE_PATH=$OPTARG ;;
        *) usage ;;
    esac
done

if [ "$EXP_PATH" == "" -o "$SUR_TIME_PATH" == "" -o "$CANCER" == "" ];
then
    usage
fi

mkdir $BASE_PATH/$CANCER
mkdir $BASE_PATH/$CANCER/PDF
mkdir $BASE_PATH/$CANCER/PDF/poor_prognosis
mkdir $BASE_PATH/$CANCER/PDF/good_prognosis
Rscript $SCRIPTS_HOME/survival.R $EXP_PATH $SUR_TIME_PATH $CANCER $LABEL $SUR_PVALUE $BOUND $BASE_PATH

