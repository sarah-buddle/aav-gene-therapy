for vector in vector pSMN pAAV29 combined; do

  # Split each read into individual fastq files
  mkdir ${results}/${sample}/${vector}_reads

  split -l 4 --numeric-suffixes --suffix-length=4 \
  ${results}/${sample}/${sample}_${vector}_mapped_filtered.fastq \
  ${results}/${sample}/${vector}_reads/${vector}_read_

  # Convert to fasta
  for file in ${results}/${sample}/${vector}_reads/*; do

      cat ${file} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > ${file}.fasta

  done

  # Move fastqs to a separate file
  mkdir ${results}/${sample}/${vector}_reads/fastq
  find ${results}/${sample}/${vector}_reads/ -type f ! -iname "*.fasta" -maxdepth 1 -exec mv -t ${results}/${sample}/${vector}_reads/fastq {} \;

  for file in ${results}/${sample}/${vector}_reads/${vector}_read*.fasta; do

      grep -v '>' $file | wc -m >> ${results}/${sample}/${vector}_reads/length.txt

  done

# Calculate lengths
  cat ${results}/${sample}/${sample}_${vector}.fastq |
  awk  'NR%4==2 {print length}' | 
  sort -nr \
  > ${results}/${sample}/${vector}_reads/length_sorted.txt 

  # Calculate total bases and N50
  sort -nr ${results}/${sample}/${vector}_reads/length.txt > ${results}/${sample}/${vector}_reads/length_sorted.txt

  TOTAL_LENGTH=$(awk '{SUM+=$1}END{print SUM}' ${results}/${sample}/${vector}_reads/length_sorted.txt)

  echo Total bases $TOTAL_LENGTH >> ${results}/${sample}/${vector}_reads/stats.txt

  CUMULATIVE_LENGTH=0

  # Initialize a variable to store the n50 value
  N50=0

  # Read through the sorted contig lengths one by one
  while read LENGTH; do
    # Increment the cumulative length by the length of the current contig
    CUMULATIVE_LENGTH=$((CUMULATIVE_LENGTH + LENGTH))

    # If the cumulative length is greater than or equal to half of the total length, then we have found the n50 value
    if [ $CUMULATIVE_LENGTH -ge $((TOTAL_LENGTH / 2)) ]; then
      N50=$LENGTH
      break
    fi
  done < ${results}/${sample}/${vector}_reads/length_sorted.txt

  echo N50 $N50 >> ${results}/${sample}/${vector}_reads/stats.txt

  # Make alignment dot plots
  for file in *.fasta; do

    redotable --window 20 ${file} \
    ${vectors}/${vector}.fasta \
    plots_vector/${file}_redotable.png

  done

done