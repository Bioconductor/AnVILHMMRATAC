task hmmratac_run {
    File bam_file
    File fastq1

    String out = basename(fastq1, "_1.fastq.gz")

    command <<<
        set -e -o pipefail
        samtools sort ${bam_file} -o ATACseq.sorted.bam
        samtools index -@ 6 ATACseq.sorted.bam ATACseq.sorted.bam.bai
        samtools view -H ATACseq.sorted.bam | \
        awk -F'[\t:]' '$1 == "@SQ" {print $3"\t"$5}' > genome.info
        java -jar /HMMRATAC_V1.2.10_exe.jar --window 250000 \
        -b ATACseq.sorted.bam -i ATACseq.sorted.bam.bai -g genome.info \
        -o ${out}
        awk -v OFS="\t" '$13 >= 10 {print}' ${out}_peaks.gappedPeak > ${out}.filteredPeaks.gappedPeak
    >>>

    runtime {
        docker: "mtmorgan/hmmratac:latest"
        memory: "10G"
        cpu: "6"
        disk: "local-disk 100 SSD"
    }

     output {
        File filtered_gappedPeak = "${out}.filteredPeaks.gappedPeak"
    }
}

workflow hmmratac {
    File bam_file
    File fastq1

    call hmmratac_run {
        input:
        bam_file = bam_file,
        fastq1 = fastq1
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
