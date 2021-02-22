task bam_align {
    File ref
    File fastq1
    File fastq2

    command {
        bwa index ${ref}
        bwa mem ${ref} ${fastq1} ${fastq2}
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }

    output {
        File ATACseq = "ATACseq.bam"
    }
}

task hmmratac_index {
    File bam_file
    String index_name

    command {
        samtools sort ${bam_file} -o ATACseq.sorted.bam
        samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai
        samtools view -H ATACseq.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > ${index_name}.info
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }

    output {
        File bam_index = "${index_name}.info"
        File sorted_bam = "ATACseq.sorted.bam"
        File sorted_bam_index = "ATACseq.sorted.bam.bai"
    }
}

task hmmratac_run {
    File sorted_bam
    File sorted_bam_index
    File bam_index
    String filter

    command {
    java -jar HMMRATAC_V1.2.4_exe.jar -b ${sorted_bam} -i ${sorted_bam_index} -g ${bam_index}
    awk -v OFS="\t" '$13>=${filter} {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }

    output {
        File hmmratac_output = NAME.filteredPead.gappedPeak
    }
}

workflow hmmratac {
    File ref
    File fastq1
    File fastq2
    String index_name
    String filter

    call bam_align {
        input:
        ref = ref,
        fastq1 = fastq1,
        fastq2 = fastq2
    }

    call hmmratac_index {
        input:
        bam_file = bam_align.ATACseq,
        index_name = index_name
    }

    call hmmratac_run {
        input:
        sorted_bam = hmmratac_index.sorted_bam,
        sorted_bam_index = hmmratac_index.sorted_bam_index,
        bam_index = hmmratac_index.bam_index,
        filter = filter
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
