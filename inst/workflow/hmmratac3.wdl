task hmmratac_run {
    File bam_file

    command <<<
        set -e -o pipefail
        samtools sort ${bam_file} -o ATACseq.sorted.bam
        samtools index -@ 6 ATACseq.sorted.bam ATACseq.sorted.bam.bai
        samtools view -H ATACseq.sorted.bam | \
        awk -F'[\t:]' '$1 == "@SQ" {print $3"\t"$5}' > genome.info
        java -jar /HMMRATAC_V1.2.10_exe.jar -b ATACseq.sorted.bam \
        -i ATACseq.sorted.bam.bai -g genome.info
        awk -v OFS="\t" '$13 >= 10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak
    >>>

    runtime {
        docker: "mtmorgan/hmmratac:latest"
        memory: "10G"
        cpu: "6"
        disk: "local-disk 100 SSD"
    }

     output {
        File filtered_gappedPeak = "NAME.filteredPeaks.gappedPeak"
    }
}

workflow hmmratac {
    File bam_file

    call hmmratac_run {
        input:
        bam_file = bam_file,
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
