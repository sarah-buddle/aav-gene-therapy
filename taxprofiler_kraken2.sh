### Running kraken2 + bracken though nf-core pipeline taxprofiler

# Illumina
nextflow run nf-core/taxprofiler -r dev -profile singularity -w ${results}/work --max_cpus 4 \
--input ${results}/samplesheets/samplesheet_${run}_test1.csv --outdir ${results} \
--databases databases.csv \
--perform_shortread_qc \
--perform_shortread_hostremoval \
--hostremoval_reference ${human_genome} \
--shortread_hostremoval_index ${human_index_dir_bowtie} \
--save_hostremoval_unmapped \
--run_kraken2 --kraken2_save_readclassifications \
--run_bracken

# ONT
 nextflow run nf-core/taxprofiler -r dev -profile conda -with-tower -w ${results}/work --max_cpus 4 \
-c ${software}/nf-core-configs/custom_resources.conf \
--input ${results}/samplesheets/samplesheet_${run}_test1.csv \
--databases ${results}/databases/database_kraken_refseq_nucleotide_v2.csv \
--outdir ${results} \
--perform_longread_qc \
--longread_qc_skipqualityfilter \
--perform_longread_hostremoval \
--hostremoval_reference ${human_genome} \
--save_hostremoval_index \
--save_hostremoval_unmapped 

# Bracken run separately
bracken -d $db \
-i ${results}/kraken2/${db_name}/${sample}_se_${run}_${db_name}.kraken2.kraken2.report.txt \
-o ${results}/bracken/${db_name}/${sample}_se_${run}_${db_name}.bracken.tsv



