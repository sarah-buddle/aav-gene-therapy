### Align trimmed reads to manufacturing plasmid sequences ###

for vector in pSMN pAAV29 pHelper; do

    # Create minimap2 index
    ${minimap2} -d ${vectors}/${vector}.mmi \
    ${vectors}/${vector}.fasta

    # Map reads using minimap2 and extract aligned reads
    ${minimap2} -ax map-ont -t 4 ${vectors}/${vector}.mmi \
    ${results}/${sample}/${sample}_porechop_alladapters.fastq |
    ${samtools} view -bF 4 -h - |
    ${samtools} sort - \
    > ${results}/${sample}/${sample}_${vector}_mapped.bam

    # Create index
    ${samtools} index ${results}/${sample}/${sample}_${vector}_mapped.bam

    # Convert to fastq
    ${samtools} fastq ${results}/${sample}/${sample}_${vector}_mapped.bam \
    > ${results}/${sample}/${sample}_${vector}.fastq

done

# Run filter_sam.R

for vector in pSMN pAAV29 pHelper; do

    $samtools view -bh -N ${results}/${sample}/${sample}_${vector}_mapped_filtered.txt \
    ${results}/${sample}/${sample}_${vector}_mapped.bam \
    > ${results}/${sample}/${sample}_${vector}_mapped_filtered.bam

    $samtools fastq ${results}/${sample}/${sample}_${vector}_mapped_filtered.bam \
    > ${results}/${sample}/${sample}_${vector}_mapped_filtered.fastq

done

