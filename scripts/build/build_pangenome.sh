#!/usr/bin/env bash

shopt -s expand_aliases
source ~/.bash_aliases

# set -e
# set -u
set -o pipefail

# source scripts/utils/setup_analysis.sh 'pectobacterium.v2'
source scripts/utils/setup_analysis.sh $@

if [ -z ${pan_db+x} ]; then
    echo "\$pan_db is unset"
    error_exit 1
fi

source $TOOLS_PATH/miniconda3/etc/profile.d/conda.sh
conda activate pantools_v4_3
export PANTOOLS="pantools -Xms40g -Xmx120g"
# conda activate pantools_master
# export PANTOOLS="$PANTOOLS_4_1"
######################################################################

#TMPDIR="$LOCAL_DIR_PATH/tmp"
#[ ! -d ${TMPDIR} ] && mkdir -p ${TMPDIR}
#[ -d ${TMPDIR}/pantools ] && rm -rd ${TMPDIR}/pantools
#[ ! -d ${TMPDIR}/pantools ] && mkdir -p ${TMPDIR}/pantools

#process_start build_pangenome
#$PANTOOLS build_pangenome --threads 40 --scratch-directory ${TMPDIR}/pantools \
#    --num-db-writer-threads 10 --cache-size 25000000 \
#    ${pan_db} $PANGENOME_DIR/genomes_fa.list
#error_exit $?

#cp -r --preserve ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.raw

######################################################################

### add annotations
#process_start add_annotations
#$PANTOOLS add_annotations --connect ${pan_db} $PANGENOME_DIR/genomes_gff3.list
#error_exit $?
#export PANTOOLS="$PANTOOLS_4_1"

### add phenotypes
#process_start add_phenotypes
#printf "y\n" | $PANTOOLS remove_phenotypes ${pan_db}
#$PANTOOLS remove_phenotypes ${pan_db}
#$PANTOOLS add_phenotypes ${pan_db} $PANGENOME_DIR/genomes_metadata.csv
#error_exit $?

### add_functions
#process_start add_InterProScan_annotations
#$PANTOOLS remove_functions ${pan_db}
#$PANTOOLS add_functions ${pan_db} $PANGENOME_DIR/functional_annotations.txt
#error_exit $?

### COG
#process_start add_COG_annotations
#$PANTOOLS_DEV add_functions ${pan_db} $PANGENOME_DIR/eggnog_annotations.txt
#error_exit $?

#### BUSCO
#process_start busco_protein
#$PANTOOLS busco_protein -t 30 --busco10 enterobacterales_odb10 ${pan_db}
#error_exit $?

#cp -r --preserve ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.fn

######################################################################
## type strains
#bash ./scripts/b_construction/grouping_subsets_process.sh "typeStrain" "374,96,385,386,379,388,116,375,377,347,250,335,338,265,256,32,269,387,266"

### generate random subsets to run `optimal_grouping`
#Rscript scripts/c_analysis/grouping_subsets_build.R

### run `optimal_grouping` on these random subsets
#nohup \
## sed -n '1,6!p' analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
#cat analysis/03_pangenome_pecto_v2/subset_optimal_group/subsets_conf.tab | \
#parallel --colsep '\t' --jobs 10  --keep-order --workdir $PWD --halt soon,fail=1 \
#--load 100% --results logs/pantools/sub_opt_group/{1} \
#--joblog logs/pantools/sub_opt_group/sub_opt_group.log \
#bash ./scripts/b_construction/grouping_subsets_process.sh {1} {2} \
#>>logs/pantools/sub_opt_group/nohup.out 2>&1 &

### summarize the results
#Rscript scripts/c_analysis/grouping_subsets_analyze.R

######################################################################
### grouping with relaxation 4 setting
#process_start group_v4
#nice $PANTOOLS group -t 30 --relaxation 4 ${pan_db}
#error_exit $?

#cp -r --preserve ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.grp

## optimized grmultiqc_buscoouping
#$PANTOOLS remove_grouping -v all ${pan_db}
#process_start optimal_grouping
#nice $PANTOOLS optimal_grouping -t 40 ${pan_db} ${pan_db}/busco/enterobacterales_odb10
#error_exit $?

#Rscript ${pan_db}/optimal_grouping/optimal_grouping.R
#$PANTOOLS grouping_overview ${pan_db}

#$PANTOOLS change_grouping -v 4 ${pan_db}
#$PANTOOLS grouping_overview ${pan_db}

#cp -rp ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.opt_grp

######################################################################
### extract the metrics from pangenome
#process_start extract_pangenome_metrics
#$PANTOOLS metrics ${pan_db}
#error_exit $?

### Core unique thresholds
#process_start core_unique_thresholds
#$PANTOOLS core_unique_thresholds ${pan_db}
#error_exit $?
#Rscript ${pan_db}/core_unique_thresholds/core_unique_thresholds.R

### gene classification: core and unique
#process_start gene_classification_core_unique
#$PANTOOLS gene_classification ${pan_db}
#error_exit $?

### Gene distance tree
#Rscript ${pan_db}/gene_classification/genome_gene_distance_tree.R

#mv ${pan_db}/gene_classification ${pan_db}/gene_classification.100.0

### gene classification: soft core and cloud
#process_start gene_classification_softcore_cloud
#$PANTOOLS gene_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
#error_exit $?

### Gene distance tree
#Rscript ${pan_db}/gene_classification/genome_gene_distance_tree.R
#mv ${pan_db}/gene_classification ${pan_db}/gene_classification.95.5

### homology group information
#sed -i.bak -r -n '/^#/! p' ${pan_db}/gene_classification.100.0/group_identifiers/all_homology_groups.csv
#process_start "extracting homology group information"
#$PANTOOLS group_info -H ${pan_db}/gene_classification.100.0/group_identifiers/all_homology_groups.csv ${pan_db}
#error_exit $?

##########################################################################
### k-mer classification
## K-mer classification: soft core and cloud
#process_start kmer_classification_core_unique
#$PANTOOLS k_mer_classification ${pan_db}
#error_exit $?

#Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
#mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.100.0

### K-mer classification: soft core and cloud
#process_start kmer_classification_softcore_cloud
#$PANTOOLS k_mer_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
#error_exit $?

#Rscript ${pan_db}/kmer_classification/genome_kmer_distance_tree.R
#mv ${pan_db}/kmer_classification ${pan_db}/kmer_classification.95.5

########################################################################
##### Pangenome structure
### Pangenome structure: genes
#process_start pangenome_structure_gene
#$PANTOOLS pangenome_structure -t 20 ${pan_db}
#error_exit $?

#Rscript ${pan_db}/pangenome_size/gene/pangenome_growth.R
#Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_or_average.R
#Rscript ${pan_db}/pangenome_size/gene/gains_losses_median_and_average.R
##Rscript ${pan_db}/pangenome_size/gene/heaps_law.R

#mv ${pan_db}/pangenome_size/gene ${pan_db}/pangenome_size/gene.pangenome

### Pangenome structure: kmer
#process_start pangenome_structure_kmer
#$PANTOOLS pangenome_structure -t 20 --kmer ${pan_db}
#error_exit $?

#Rscript ${pan_db}/pangenome_size/kmer/pangenome_growth.R

######################################################################
### functional_classification: core and unique
#process_start functional_classification_core_unique
#$PANTOOLS functional_classification ${pan_db}
#error_exit $?

#mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.100.0

### functional_classification: softcore and cloud
#process_start functional_classification_softcore_cloud
#$PANTOOLS functional_classification --core-threshold 95 --unique-threshold 5 ${pan_db}
#error_exit $?

#mv ${pan_db}/function/functional_classification ${pan_db}/function/functional_classification.95.5

### function_overview
#process_start function_overview
#$PANTOOLS function_overview ${pan_db}
#error_exit $?

#Rscript ${pan_db}/function/cog_per_class.R

### GO enrichment for core, accessory and unique homology groups
#for grp in core accessory unique
#do
#    process_start GO_enrichment:$grp
#    $PANTOOLS go_enrichment -H ${pan_db}/gene_classification.100.0/${grp}_groups.csv  ${pan_db}
#    error_exit $?
#    mv ${pan_db}/function/go_enrichment ${pan_db}/function/go_enrichment.100.0.${grp}
#done

### GO enrichment for soft-core, accessory and cloud homology groups
#for grp in core accessory unique
#do
#    process_start GO_enrichment:$grp
#    $PANTOOLS go_enrichment -H ${pan_db}/gene_classification.95.5/${grp}_groups.csv  ${pan_db}
#    error_exit $?
#    mv ${pan_db}/function/go_enrichment ${pan_db}/function/go_enrichment.95.5.${grp}
#done

#cp -r --preserve ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.chr
#######################################################################

## MSA for homology groups
process_start "msa for homology groups"
$PANTOOLS msa -t 40 --method per-group --mode nucleotide ${pan_db}
error_exit $?

cp -r --preserve ${pan_db} $PANGENOME_DIR/backup/${PANGENOME_NAME}.DB.msa
#######################################################################

### SNP tree using core gene SNPs
#process_start core_phylogeny
#cp -r --preserve ${pan_db}/gene_classification.100.0 ${pan_db}/gene_classification
#$PANTOOLS core_phylogeny -t 30  --clustering-mode ML ${pan_db}
#error_exit $?
#rm -r ${pan_db}/gene_classification

######################################################################
