task hmmratac_run {
    File bwa_ref
    File bwa_ref_amb
    File bwa_ref_ann
    File bwa_ref_bwt
    File bwa_ref_pac
    File bwa_ref_sa
    File fastq1
    File fastq2

    command {
        bwa mem ${bwa_ref} ${fastq1} ${fastq2}
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }
}

workflow hmmratac {
    File bwa_ref
    File bwa_ref_amb
    File bwa_ref_ann
    File bwa_ref_bwt
    File bwa_ref_pac
    File bwa_ref_sa
    File fastq1
    File fastq2

    call hmmratac_run {
        input:
        bwa_ref = bwa_ref,
        bwa_ref_amb = bwa_ref_amb,
        bwa_ref_ann = bwa_ref_ann,
        bwa_ref_bwt = bwa_ref_bwt,
        bwa_ref_pac = bwa_ref_pac,
        bwa_ref_sa = bwa_ref_sa,
        fastq1 = fastq1,
        fastq2 = fastq2
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
} 
