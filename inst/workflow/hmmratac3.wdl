task hmmratac_run {
    File bam_file
    File genome_info
    String filter

    command {
        samtools sort ${bam_file} -o ATACseq.sorted.bam
        samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai
        java -jar HMMRATAC_V1.2.4_exe.jar -b ATACseq.sorted.bam \
        -i ATACseq.sorted.bam.bai -g ${genome_info}
        awk -v OFS="\t" '$13>=${filter} {print}' \
        NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }

     output {
        File filtered_gappedPeak = "NAME.filteredPeaks.gappedPeak"
    }
}

workflow hmmratac {
    File bam_file
    File genome_info
    String filter

    call hmmratac_run {
        input:
        bam_file = bam_file,
        genome_info = genome_info,
        filter = filter
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
