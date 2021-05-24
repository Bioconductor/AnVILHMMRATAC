task bam_index {
    File ref

    command <<<
        set -e -o pipefail
        bwa index ${ref}
        #mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        #"SELECT chrom, size FROM hg19.chromInfo" > genome_info
    >>>

    runtime {
        docker: "mtmorgan/hmmratac:latest"
        memory: "10G"
        cpu: "6"
    }

    output {
        File bwa_ref_amb = "${ref}.amb"
        File bwa_ref_ann = "${ref}.ann"
        File bwa_ref_bwt = "${ref}.bwt"
        File bwa_ref_pac = "${ref}.pac"
        File bwa_ref_sa = "${ref}.sa"
        #File genome_info = "genome_info"
    }
}

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

    String out = basename(fastq1, "_1.fastq.gz")    

    command <<<
        set -e -o pipefail
        bwa mem -t 6 ${bwa_ref} ${fastq1} ${fastq2} | \
        samtools view -bS -t ${chromInfo} -o bam_file.bam - 
        samtools sort bam_file.bam -o ${out}.sorted.bam
        samtools index -@ 6 ${out}.sorted.bam ${out}.sorted.bam.bai
        samtools view -H ${out}.sorted.bam | \
        awk -F'[\t:]' '$1 == "@SQ" {print $3"\t"$5}' > genome.info
        java -jar /HMMRATAC_V1.2.10_exe.jar --window 25000 \
        -b ${out}.sorted.bam -i ${out}.sorted.bam.bai -g genome.info \
        -o ${out}
        awk -v OFS="\t" '$13 >= 10 {print}' ${out}_peaks.gappedPeak > ${out}.filteredPeaks.gappedPeak
    >>>

    runtime {
        docker: "mtmorgan/hmmratac:latest"
        memory: "40G"
        cpu: "6"
        disks: "local-disk 100 SSD"
    }

    output {
        File filtered_gappedPeak = "${out}.filteredPeaks.gappedPeak"
        File sorted_bam = "${out}.sorted.bam"
        File sorted_bam_bai = "${out}.sorted.bam.bai"
        File summits_bed = "${out}_summits.bed"
    }
}

workflow hmmratac {
    File ref
    Array[File] fastq1
    Array[File] fastq2
    File chromInfo

    Array[Pair[File,File]] fastq_pairs = zip(fastq1, fastq2)

    call bam_index {
        input:
        ref = ref
    }

    scatter (fastq_pair in fastq_pairs) {
        File fastq_1 = fastq_pair.left
        File fastq_2 = fastq_pair.right

        File bwa_ref_amb = bam_index.bwa_ref_amb
        File bwa_ref_ann = bam_index.bwa_ref_ann
        File bwa_ref_bwt = bam_index.bwa_ref_bwt
        File bwa_ref_pac = bam_index.bwa_ref_pac
        File bwa_ref_sa = bam_index.bwa_ref_sa
        #File genome_info = bam_index.genome_info

        call hmmratac_run {
            input:
            bwa_ref = ref,
            bwa_ref_amb = bwa_ref_amb,
            bwa_ref_ann = bwa_ref_ann,
            bwa_ref_bwt = bwa_ref_bwt,
            bwa_ref_pac = bwa_ref_pac,
            bwa_ref_sa = bwa_ref_sa,
            chromInfo = chromInfo,
            fastq1 = fastq_1,
            fastq2 = fastq_2
        }
    }

    meta {
        author: "Kayla Interdonato"
        email: "Kayla.Morrell@roswellpark.org"
        description: "Running the HMMRATAC workflow."
    }
}
