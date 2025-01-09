### Align taxprofiler filtered reads to viral and manufacturing plasmid sequences ###

for genome in NC_001401 NC_000898.1 NC_001405 pSMN pAAV29 pHelper; do

    # Build index
    ${bowtie}-build ${genome_dir}/${genome}.fasta ${genome_dir}/${genome}

    # Align with bowtie2
    (${bowtie} -x ${genome_dir}/${genome}  \
    --very-sensitive \
    ##--score-min L,0,-0.1 -N 0 -L 22 --mp 6,2 --rdg 5,3 --rfg 5,3 \ # used for AAV2 WT genome
    -1 ${filtered_dir}/${sample}_${dnarna}_filtered_1.fastq \
    -2 ${filtered_dir}/${sample}_${dnarna}_filtered_2.fastq \
    -S ${results}/${dnarna}/${sample}_${dnarna}_${genome}_v2.sam) \
    2> ${results}/${dnarna}/${sample}_${dnarna}_${genome}_log_v2.txt

    # Extract aligned reads
    ${samtools} view -bF 4 -h ${results}/${dnarna}/${sample}_${dnarna}_${genome}_v2.sam |
    ${samtools} sort \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_v2.bam

    rm ${results}/${dnarna}/${sample}_${dnarna}_${genome}_v2.sam

    # Remove PCR duplicates
    $samtools collate -@ 4 -O -u ${results}/${dnarna}/${sample}_${dnarna}_${genome}_v2.bam |
    $samtools fixmate -@ 4 -m -u - - |
    $samtools sort -@ 4 -u - |
    $samtools markdup -@ 4 -r -d 2500 - ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.bam \
    -f ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.txt

    # Index
    $samtools index ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.bam

    # Calculate depth
    $samtools depth ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.bam \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_depth_v2.txt

    # Create fastq
    $samtools fastq ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.bam \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_v2.fastq

    # Create alignment plot
    $rscript $alignment_r \
    ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_depth_v2.txt \
    ${results}/${dnarna}/${sample}_${dnarna}_${genome}_dedup_depth_v2.png

    # Repeat above with raw reads
    (${bowtie} -x ${genome_dir}/${genome} --very-sensitive \
    -1 ${raw_dir}/${sample}_${dnarna}_1.fq.gz \
    -2 ${raw_dir}/${sample}_${dnarna}_2.fq.gz \
    -S ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw.sam) \
    2> ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_log.txt

    ${samtools} view -bF 4 -h ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw.sam |
    ${samtools} sort \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw.bam

    rm ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw.sam

    $samtools collate -@ 4 -O -u ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw.bam |
    $samtools fixmate -@ 4 -m -u - - |
    $samtools sort -@ 4 -u - |
    $samtools markdup -@ 4 -r -d 2500 - ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.bam \
    -f ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.txt

    $samtools index ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.bam

    $samtools depth ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.bam \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup_depth.txt

    $samtools fastq ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.bam \
    > ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup.fastq

    $rscript $alignment_r \
    ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup_depth.txt \
    ${results}/${dnarna}/${sample}_${dnarna}_${genome}_raw_dedup_depth.png

done



