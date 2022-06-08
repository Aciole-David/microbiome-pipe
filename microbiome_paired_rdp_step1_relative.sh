#!/bin/bash
#SBATCH -n 24                #Neste caso estão sendo pedidos 64 cores, configure como desejar. 
#SBATCH --ntasks-per-node=24 #Para jobs de até 64 cores manter o mesmo número que o numero de cores pedidos, para jobs com mais de 64 cores, mantter sempre em 64.
#SBATCH -p fast           #Pedido para fila média, se desejar outra fila, modifique
#SBATCH --time=24:00:00
#SBATCH --mem=64000
#SBATCH --job-name=microbiome   #Mudar o JOBNAME para o nome desejado para o job.



###############################################################################
# Experiment details
###############################################################################

###############################################################################
# Environment variables
###############################################################################
#WORKDIR=$HOME/Metagenomics/workdir
#MAP=$WORKDIR/map/microbiome_TB10.txt
#MICROBIOME=$WORKDIR/microbiome_paired
#MICOUT=$MICROBIOME/output/rdp
#MICONFIG=$WORKDIR/config/qiime_16s_parameters.txt
#MISAMPLEP=$WORKDIR/Cachexia/16S/samples
#GG=$HOME/.bin/qiime_dbs/greengenes/gg_13_8_otus/rep_set/97_otus.fasta
#UNITY=$HOME/.bin/qiime_dbs/unite/2/sh_refs_qiime_ver7_97_01.12.2017.fasta

WORKDIR=$HOME/Metagenomics/workdir/
MAP=$WORKDIR/Caquexia/map/microbiome_TB10.txt
MICROBIOME=$WORKDIR/microbiome_paired
MICOUT=$MICROBIOME/output/rdp/TB-8/
MICONFIG=$HOME/Metagenomics/config/qiime_parameters_16S.txt 
MISAMPLEP=$WORKDIR/Caquexia/16S/samples
GG=$HOME/.bin/qiime_dbs/greengenes/gg_13_8_otus/rep_set/97_otus.fasta
UNITY=$HOME/.bin/qiime_dbs/unite/2/sh_refs_qiime_ver7_97_01.12.2017.fasta


source activate qiime

echo "###############################"
echo "###############################"
echo "INICIANDO PRE-PROCESSAMENTO..."
echo "###############################"
echo "###############################"


###############################################################################
# Trimming
###############################################################################

echo " "
echo "AVISO: INICIANDO O PROCESSO DE TRIMMING DAS BIBLIOTECAS..."

mkdir -p $MICOUT/trimming
cd $MISAMPLEP

for i in $(ls | rev | cut -c 11- | rev | uniq)

do

seqtk trimfq -b 1 -e 30 ${i}_001.fastq > $MICOUT/trimming/${i}_001.fastq

done

echo "PROCESSO DE TRIMMING FINALIZADO."
echo " "

###############################################################################

###############################################################################
# Join Paired-End 
###############################################################################

echo "AVISO: INICIANDO O PROCESSO DE JOIN PAIRED-END..."

multiple_join_paired_ends.py \
                -i $MICOUT/trimming/ \
                -o $MICOUT/joined/ \
                -p $MICONFIG

echo "PROCESSO DE JOINED PAIRED-END FINALIZADO."
echo " "

###############################################################################

###############################################################################
# Renomeando pastas
###############################################################################

cd $MICOUT/joined/

for sample_paired in $(ls)

do

sample_rename=`echo ${sample_paired} | cut -d "_" -f1`
mv ${sample_paired} ${sample_rename}

done

###############################################################################

###############################################################################
# Removendo arquivos 
###############################################################################

rm */fastqjoin.un*.fastq

###############################################################################

###############################################################################
# Removendo bibliotecas zeradas
###############################################################################

mkdir $MICOUT/null
cd $MICOUT/joined/
zero=`find -size "0k" -exec ls {} \; | cut -d "/" -f2`
mv $zero $MICOUT/null

###############################################################################

###############################################################################

###############################################################################
# Demultiplex .fastq sequence data. 
# "Turning off" filter parameters, and storing the demultiplexed .fastq file. 
###############################################################################

echo "AVISO: INICIANDO O PROCESSO DE SPLITTING..." 

multiple_split_libraries_fastq.py \
                -i $MICOUT/joined/ \
                -o $MICOUT/splitting/ \
                -p $MICONFIG \
                --remove_filepath_in_name \
                --include_input_dir_path \
                --demultiplexing_method sampleid_by_file

echo "PROCESSO DE SPLITTING FINALIZADO."
echo " "

###############################################################################

###############################################################################
# Split FASTA file
###############################################################################

echo "AVISO: DIVIDINDO ARQUIVO FASTA EM 3 PARTES PARA REDUZIR O TAMANHO E PERMITIR A EXECUÇÃO PELO USEARCH61..."  


mkdir -p $MICOUT/split_fasta
cd $MICOUT/split_fasta

#pyfasta split -n3 $MICOUT/splitting/seqs.fna
#cp $MICOUT/fastx_filter/reads.*.fa .

fasta-splitter --n-parts 3 $MICOUT/splitting/seqs.fna --out-dir $MICOUT/split_fasta

echo "FINALIZADA DIVISÃO."

echo " "

###############################################################################

###############################################################################
# Identify and filter chimeric seqs
###############################################################################

echo "AVISO: IDENTIFICANDO E REMOVENDO SEQUÊNCIAS QUIMÉRICAS..."

for sample_paired in $(ls)

do

mkdir -p $MICOUT/no_chimera/$sample_paired

identify_chimeric_seqs.py \
            -m usearch61 \
            -i $MICOUT/split_fasta/$sample_paired \
            -r $GG \
            -o $MICOUT/no_chimera/$sample_paired/

filter_fasta.py \
            -f $MICOUT/split_fasta/$sample_paired \
            -o $MICOUT/no_chimera/$sample_paired/seqs_chimeras_filtered.fna \
            -s $MICOUT/no_chimera/$sample_paired/chimeras.txt -n

echo "SEQUÊNCIAS QUIMÉRICAS REMOVIDAS."
echo " "

done

###############################################################################

###############################################################################
# Joined FASTA file
###############################################################################

mkdir -p $MICOUT/joined_fasta
cd $MICOUT/joined_fasta

cat $MICOUT/no_chimera/*/seqs_chimeras_filtered.fna > $MICOUT/joined_fasta/seqs_chimeras_filtered.fna

###############################################################################

echo "###############################"
echo "###############################"
echo "PRE-PROCESSAMENTO FINALIZADO..."
echo "###############################"
echo "###############################"

echo " "

echo "###############################"
echo "###############################"
echo "INICIANDO PROCESSAMENTO..."
echo "###############################"
echo "###############################"

###############################################################################
# PICKING OTU
###############################################################################

echo " "
echo "AVISO: INICIANDO O PROCESSO DE PICK OTU USANDO SORTMERNA..."

mkdir -p $MICOUT/otu_picking/

pick_open_reference_otus.py \
                -m sortmerna_sumaclust \
                -i $MICOUT/joined_fasta/seqs_chimeras_filtered.fna \
                -o $MICOUT/otu_picking/rdp \
                -r $GG \
                -p $MICONFIG \
                --force
         
cp $MICOUT/otu_picking/rdp/otu_table_mc2_w_tax_no_pynast_failures.biom \
         $MICOUT/otu_picking/rdp/otu_table_sortmerna.biom

echo "PICK OTU FINALIZADO."
echo " "
echo "AVISO: AJUSTANDO OTU TABLE..."


filter_taxa_from_otu_table.py \
		-i $MICOUT/otu_picking/rdp/otu_table_sortmerna.biom \
		-o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
		-n c__Chloroplast,f__mitochondria

filter_otus_from_otu_table.py \
		-i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
		-o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
		-n 2

## Abundance 0001


filter_otus_from_otu_table.py \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                --min_count_fraction 0.0001


filter_fasta.py \
                -f $MICOUT/otu_picking/rdp/rep_set.fna \
                -b $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MICOUT/otu_picking/rdp/rep_set_0001.fna

assign_taxonomy.py \
                -o $MICOUT/otu_picking/rdp/rdp_assigned_taxonomy_0001 \
                -i $MICOUT/otu_picking/rdp/rep_set_0001.fna \
                -p $MICONFIG
                
align_seqs.py \
                -i $MICOUT/otu_picking/rdp/rep_set_0001.fna \
                -o $MICOUT/otu_picking/rdp/pynast_aligned_seqs_0001/ \
                -t $HOME/.scripts/qiime_dbs/greengenes/gg_13_8_otus/rep_set_aligned/97_otus.fasta \
                --alignment_method pynast \
                --pairwise_alignment_method uclust \
                --min_percent_id 0.75

filter_alignment.py \
                -o $MICOUT/otu_picking/rdp/pynast_aligned_seqs_0001 \
                -i $MICOUT/otu_picking/rdp/pynast_aligned_seqs_0001/rep_set_0001_aligned.fasta
                
make_phylogeny.py \
                -i $MICOUT/otu_picking/rdp/pynast_aligned_seqs_0001/rep_set_0001_aligned_pfiltered.fasta \
                -o $MICOUT/otu_picking/rdp/rep_set_0001.tre \
                --root_method tree_method_default \
                --tree_method fasttree

## Abundance 00001

filter_otus_from_otu_table.py \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                --min_count_fraction 0.00001


filter_fasta.py \
                -f $MICOUT/otu_picking/rdp/rep_set.fna \
                -b $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MICOUT/otu_picking/rdp/rep_set_00001.fna

assign_taxonomy.py \
                -o $MICOUT/otu_picking/rdp/rdp_assigned_taxonomy_00001 \
                -i $MICOUT/otu_picking/rdp/rep_set_00001.fna \
                -p $MICONFIG
                
align_seqs.py \
                -i $MICOUT/otu_picking/rdp/rep_set_00001.fna \
                -o $MICOUT/otu_picking/rdp/pynast_aligned_seqs_00001/ \
                -t $HOME/.scripts/qiime_dbs/greengenes/gg_13_8_otus/rep_set_aligned/97_otus.fasta \
                --alignment_method pynast \
                --pairwise_alignment_method uclust \
                --min_percent_id 0.75

filter_alignment.py \
                -o $MICOUT/otu_picking/rdp/pynast_aligned_seqs_00001 \
                -i $MICOUT/otu_picking/rdp/pynast_aligned_seqs_00001/rep_set_00001_aligned.fasta
                
make_phylogeny.py \
                -i $MICOUT/otu_picking/rdp/pynast_aligned_seqs_00001/rep_set_00001_aligned_pfiltered.fasta \
                -o $MICOUT/otu_picking/rdp/rep_set_00001.tre \
                --root_method tree_method_default \
                --tree_method fasttree

echo "AJUSTES DA OTU TABLE FINALIZADOS."
echo " "

echo "AVISO: CONVERTENDO OTU TABLE PARA TSV..."

biom convert \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv

biom convert \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna.biom \
                -o $MICOUT/otu_picking/rdp/otu_table_sortmerna.txt  \
                -m $MAP  \
                --header-key=taxonomy --output-metadata-id="Consensus Lineage" \
                --process-obs-metadata=taxonomy \
                --table-type="OTU table" \
                --to-tsv


echo "CONVERSÃO DA OTU TABLE PARA TSV FINALIZADA."
echo " "

echo "AVISO: GERANDO RELATÓRIO DA OTU TABLE..."

biom summarize-table \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
                -o $MICOUT/otu_picking/rdp/biom_summary_0001.txt

biom summarize-table \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
                -o $MICOUT/otu_picking/rdp/biom_summary_00001.txt


biom summarize-table \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
                -o $MICOUT/otu_picking/rdp/biom_summary_singleton.txt

biom summarize-table \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
                -o $MICOUT/otu_picking/rdp/biom_summary_filtered.txt


biom summarize-table \
                -i $MICOUT/otu_picking/rdp/otu_table_sortmerna.biom \
                -o $MICOUT/otu_picking/rdp/biom_summary.txt


echo "RELATÓRIO FINALIZADO."
echo " "

echo "###############################"
echo "###############################"
echo "PROCESSAMENTO FINALIZADO..."
echo "###############################"
echo "###############################"


echo " "
echo " "

echo "###############################"
echo "###############################"
echo "PROCESSAMENTO FINALIZADO..."
echo "###############################"
echo "###############################"

echo " "

###############################################################################
# Core Microbiome
###############################################################################


echo "###############################"
echo "###############################"
echo "CRIANDO CORE MICROBIOME..."
echo "###############################"
echo "###############################"


echo " "
echo " "

mkdir -p $MICOUT/otu_picking/rdp/core

### 0001 abundance

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_CAQ-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_SUP-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001_TB-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_100.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_CAQ-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_SUP-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_0001_TB-B/core_table_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom


biom convert \
	-i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom \
	-o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.txt \
	-m $MAP \
	--header-key=taxonomy \
	--output-metadata-id="Consensus Lineage" \
	--process-obs-metadata=taxonomy \
	--table-type="OTU table" \
	--to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv


### 00001 abundance

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'


compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_CAQ-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_SUP-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B


compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001_TB-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_100.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_CAQ-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_SUP-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_00001_TB-B/core_table_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom


biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

### Singletons

### singleton abundance

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-B.biom \
        -m $MAP \
        -s 'ALL:CAQ-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-B.biom \
        -m $MAP \
        -s 'ALL:SUP-B'

filter_samples_from_otu_table.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -o $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-B.biom \
        -m $MAP \
        -s 'ALL:TB-B'

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_CAQ-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_SUP-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B

compute_core_microbiome.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_TB-B.biom \
        -o $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_100.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_100.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_90.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_80.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom

merge_otu_tables.py \
        -i $MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_CAQ-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_SUP-B/core_table_70.biom,$MICOUT/otu_picking/rdp/core/otu_table_core_abundance_singleton_TB-B/core_table_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom


biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv

biom convert \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
        -o $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.txt \
        -m $MAP \
        --header-key=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        --process-obs-metadata=taxonomy \
        --table-type="OTU table" \
        --to-tsv


###############################################################################
# Make OTU network
###############################################################################

echo "###############################"
echo "###############################"
echo "CRIANDO OTU NETWORK..."
echo "###############################"
echo "###############################"


echo " "
echo " "

mkdir -p $MICOUT/otu_picking/rdp/networks/core


make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/otu_network_singleton

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/otu_network_filtered

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/otu_network

### 0001 abundance

make_otu_network.py \
	-i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_0001.biom \
	-m $MAP \
	-o $MICOUT/otu_picking/rdp/networks/otu_network_0001


make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_100.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_90.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_90

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_80.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_80

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_0001_70.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_70

### 00001 abundance

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/otu_table_sortmerna_filtered_singleton_abundance_00001.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/otu_network_00001

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_100.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_90.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_90

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_80.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_80

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_abundance_00001_70.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_70

### Core Singletons

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_100.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_0001_100

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_90.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_90

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_80.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_80

make_otu_network.py \
        -i $MICOUT/otu_picking/rdp/core/merged_otu_table_core_singleton_70.biom \
        -m $MAP \
        -o $MICOUT/otu_picking/rdp/networks/core/otu_network_core_abundance_00001_70


