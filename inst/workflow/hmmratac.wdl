# need a task for BAM file aligner
# output will be ATACseq.bam

task hmmratac_index {


    command {
        samtools sort ATACseq.bam -o ATACseq.sorted.bam
        samtools index ATACseq.sorted.bam ATACseq.sorted.bam.bai
        samtools view -H ATACseq.sorted.bam| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > genome.info
    }

    runtime {
        docker:
    }

    output {
        File
    }
}

task hmmratac_run {

    command {
    java -jar HMMRATAC_V1.2.4_exe.jar -b ATACseq.sorted.bam -i ATACseq.sorted.bam.bai -g genome.info
    awk -v OFS="\t" '$13>=10 {print}' NAME_peaks.gappedPeak > NAME.filteredPeaks.gappedPeak
    }

    runtime {
        docker:
    }

    output {
        File
    }
}

workflow hmmratac {
    

    call hmmratac_index {
        input:
        
    }

    call hmmratac_run {
        input:
        
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
