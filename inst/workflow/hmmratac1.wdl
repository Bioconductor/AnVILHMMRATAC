task bam_align {
    File ref

    command {
        bwa index ${ref}
    }

    runtime {
        docker: "mtmorgan/hmmratac:latest"
    }
}

workflow hmmratac {
    File ref

    call bam_align {
        input:
        ref = ref
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
