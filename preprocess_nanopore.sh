mkdir -p ${results}/${sample}/aav2_reads

# Combine files
zcat ${data}/*.fastq.gz > ${data}/${sample}.fastq

gzip ${data}/${sample}.fastq

# Trim adaptors
mamba activate porechop

porechop -i ${data}/${sample}.fastq.gz --adapter_threshold 85 \
-o ${results}/${sample}/${sample}_porechop_alladapters.fastq -v 2 |
tee -a > ${results}/${sample}/${sample}_porechop_alladapters.log

gzip ${results}/${sample}/${sample}_porechop_alladapters.fastq

# Filter human reads
${minimap2} -ax map-ont -t 8 ${human_genome}.mmi \
${results}/${sample}/${sample}_porechop_alladapters.fastq.gz |
${samtools} view -bh - > ${results}/${sample}/${sample}_human.bam

${samtools} view -bf 4 -h ${results}/${sample}/${sample}_human.bam \
>${results}/${sample}/${sample}_human_mapped.bam

${samtools} fastq ${results}/${sample}/${sample}_human_mapped.bam \
> ${results}/${sample}/${sample}_human_filtered.fastq