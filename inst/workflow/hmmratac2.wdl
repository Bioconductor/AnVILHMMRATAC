task hmmratac_run {
    File bwa_ref
    File bwa_ref_amb
    File bwa_ref_ann
    File bwa_ref_bwt
    File bwa_ref_pac
    File bwa_ref_sa
    File chromInfo
    File fastq1
    File fastq2

    command {
        bwa mem -t 6 ${bwa_ref} ${fastq1} ${fastq2} | \
        samtools view -bS -t ${chromInfo} -o bam_file.bam - 
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
        memory: "10G"
        cpu: "6"
        disks: "local-disk 100 SSD"
    }

    output {
        File bam_file = "bam_file.bam"
    }
}

workflow hmmratac {
    File bwa_ref
    File bwa_ref_amb
    File bwa_ref_ann
    File bwa_ref_bwt
    File bwa_ref_pac
    File bwa_ref_sa
    File chromInfo
    Array[File] fastq1
    Array[File] fastq2

    scatter (fastqs in zip(fastq1, fastq2)) {
        call hmmratac_run {
            input:
            bwa_ref = bwa_ref,
            bwa_ref_amb = bwa_ref_amb,
            bwa_ref_ann = bwa_ref_ann,
            bwa_ref_bwt = bwa_ref_bwt,
            bwa_ref_pac = bwa_ref_pac,
            bwa_ref_sa = bwa_ref_sa,
            chromInfo = chromInfo,
            fastq1 = fastqs.left,
            fastq2 = fastqs.right
        }
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
} 
