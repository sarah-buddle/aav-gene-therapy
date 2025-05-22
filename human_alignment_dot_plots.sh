# Extract human mapped reads from bam file produced in preprocessing
${samtools} view -bF 4 -h ${results}/${sample}/${sample}_human_sorted.bam \
> ${results}/${sample}/${sample}_human_mapped.bam

$samtools depth -a ${results}/${sample}/${sample}_human_mapped.bam \
> ${results}/${sample}/${sample}_human_mapped_depth.txt

$samtools index ${results}/${sample}/${sample}_human_mapped.bam

# Extract ACTB from total alignment
$samtools view -bh ${results}/${sample}/${sample}_human_mapped.bam "7:5,526,409-5,563,902" \
> ${results}/${sample}/${sample}_human_mapped_beta_actin.bam

# Extract GTF2H2 gene
$samtools view -bh ${results}/${sample}/${sample}_human_mapped.bam "5:71,035,582-71,067,652" \
> ${results}/${sample}/${sample}_human_mapped_beta_gtf2h2.bam

for gene in beta_actin gtf2h2; do

    $samtools index ${results}/${sample}/${sample}_human_mapped_${gene}.bam

    $samtools fastq ${results}/${sample}/${sample}_human_mapped_${gene}.bam \
    > ${results}/${sample}/${sample}_human_mapped_${gene}.fastq

    NanoPlot -t 4 --huge --fastq ${results}/${sample}/${sample}_human_mapped_${gene}.fastq \
    -o ${results}/${sample}/${gene}_nanoplot \
    --only-report

    mkdir ${results}/${sample}/${gene}_reads

    lines=$(wc -l ${results}/${sample}/${sample}_human_mapped_${gene}.fastq | awk '{print $1}')
    max_reads=$(echo $(( $lines / 4)))

    for ((i=1; i<=max_reads; i++)); do

        let "n_reads = $i * 4"

        head -$n_reads ${results}/${sample}/${sample}_human_mapped_${gene}.fastq | tail -4 | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" \
        > ${results}/${sample}/${gene}_reads/${sample}_${gene}_read${i}.fasta

        redotable --window 20 ${results}/${sample}/${gene}_reads/${sample}_${gene}_read${i}.fasta \
        ${human_genome}/human_${gene}_region.fasta \
        ${results}/${sample}/${gene}_reads/plots/${sample}_${gene}_read${i}_redotable.png

    done

done
