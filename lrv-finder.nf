#!/usr/bin/env nextflow

/*
   ------------------------------------------------------------------------------------------------
   WELCOME TO THE lrv-finder PIPELINE
   ------------------------------------------------------------------------------------------------
   Usage:

   nextflow run lrv-finder.nf --sra "SRR123456"

   optional parameters:
   --output "output/dir/" [default="current-directory/output/"]
   --threads # (number of CPUs) [default=4]
   --memory # (amount of RAM) [in GB, default=8]
   --tempdir "tmp/dir" [default="/tmp"]
   --referenceseqs "reference_sequences.fasta" [default = "resources/lrv-genomes.fasta"]
   ------------------------------------------------------------------------------------------------
   Contact:

   Austin R. Manny
   Nibert Lab @ Harvard Medical School
   austinmanny@g.harvard.edu
   github.com/austinreidmanny
   ------------------------------------------------------------------------------------------------
   Complete background, project objectives, and information available @
   github.com/austinreidmanny/lrv-finder/README.md
   ------------------------------------------------------------------------------------------------
*/

// ---------------------------------------------------------------------------------------------- //
// Welcome the user
// ---------------------------------------------------------------------------------------------- //

println "\n====================================================================================\n" +
        "Welcome to the 'lrv-finder' pipeline! To run it, just type the following: \n\n" +
        "   nextflow run lrv-finder.nf --sra 'SRA00001,SRA00002,SRA00003' \n\n" +
        "For detailed information on this pipeline, refer to the 'README.md' file or visit \n" +
        "www.github.com/austinreidmanny/lrv-finder/ \n" +
        "====================================================================================\n"

// ---------------------------------------------------------------------------------------------- //

// ---------------------------------------------------------------------------------------------- //
// Main code
// ---------------------------------------------------------------------------------------------- //

// Create a 'run_name' parameter for naming files;
//     if single SRA accession given, name it that; if multple SRAs: "SRA_first-SRA_last"
if (params.sra.split(",").size() > 1) {

    run_name = params.sra.split(",")[0] +
               "-" +
               params.sra.split(",")[-1]
} else {
    run_name = params.sra
}

process log_inputs {
    // Save all inputs used for this run; important because files are named according to 'run_name'
    // which only indicates the first and last sample

    publishDir "${params.output}/00_analysis_info/", mode: "copy"

    input:
    val samples from params.sra

    output:
    file "${run_name}.readme.txt"

    """
    echo "Pipeline began at: \$(date)" > \
         "${run_name}.readme.txt"

    echo "Input samples: $samples" >> \
         "${run_name}.readme.txt"

    echo "Reference sequences to map to: $params.referenceseqs" >> \
         "${run_name}.readme.txt"
    """

}

process parse_sra_ids {
    // Process the SRA accessions provided by the user

    input:
    val sra_id from params.sra.split(",")

    output:
    val sra_id into sra_accessions

    """
    echo "Preparing run for SRA ID: $sra_id ..."
    """

}

process download_sra_files {
    // Take in each SRA accession number, download the files, and send them to the mapping process

    //publishDir "${params.output}/01_reads/", mode: "copy"

    input:
    val sra_id from sra_accessions

    output:
    tuple val(sra_id), file("${sra_id}*fastq") into sra_fastqs

    """
    download_sra.sh -s $sra_id -t $params.tempdir -m $task.memory -n $task.cpus -o ./
    """

}

process mapping {
    // Take in FASTQ(s), put out mapped-cleaned-sorted BAM
    publishDir "${params.output}/02_mapped_bam/", mode: "copy"

    input:
    tuple val(sra_id), file(reads) from sra_fastqs

    output:
    file "*bam" into mapped_bam
    file "*mapping-summary.txt"

    """
    mapping.sh $sra_id $params.referenceseqs $reads $task.cpus
    """

}

process combine_bams {

    publishDir "${params.output}/03_combined_bam/", mode: "copy"

    input:
    file mapped_bams from mapped_bam.collect()

    output:
    tuple val(run_name), file("${run_name}.combined.bam") into combined_bam_for_fasta
    file "${run_name}.combined.bam" into combined_bam

    """
    samtools merge ${run_name}.combined.bam $mapped_bams
    """
}

process bam_to_fasta {
    // Take in mapped BAM, split it up by viral species & save as species-specific FASTA files
    publishDir "${params.output}/04_binned_reads/", mode: "copy"

    input:
    tuple val(run_name), file(bam) from combined_bam_for_fasta

    output:
    file "*fasta.gz" into mapped_species
    file "*fasta.gz"

    """
    # The third column of BAM files is the reference it's mapped to; extract all unique references
    samtools view $bam | cut -f 3 | sort -u > "references_with_mapped_reads.txt"

    # Index the BAM so we can use the reference sequences
    samtools index $bam

    # Loop through all reference sequences with mapped reads, and extract those binned reads
    while read -r reference; do

        samtools view -b $bam \$reference | \
        samtools fasta -N - > "${run_name}.reads-mapped-to.\${reference}.fasta"
        gzip "${run_name}.reads-mapped-to.\${reference}.fasta"

        done < "references_with_mapped_reads.txt"
    """

}

process count_mapped_reads {
    // Take in the combined BAM file and retrieve the mapping counts per reference sequence

    publishDir "${params.output}/05_results_table/", mode: "copy"

    input:
    file combined_bam from combined_bam

    output:
    file "${run_name}.initial_table.txt"
    file "${run_name}.initial_table.txt" into initial_table

    """
    # Get mapping statistics, eliminate uninformative columns, then remove the last (unhelpful) line

    samtools idxstats $combined_bam |
    cut -f 1,3 |
    awk '\$1 != "*" {print \$0}' > \
    "${run_name}.initial_table.txt"
    """

}


process generate_final_summary{
    // Take in the rough table and create a final summary

    publishDir "${params.output}/05_results_table/",
        mode: "copy",
        pattern: "final_results.txt",
        saveAs: { filename -> "${run_name}.${filename}" }

    input:
    file initial_table

    output:
    file "final_results.txt" into final_table

    """
    generate_summary.sh $initial_table
    """

}

process print_results {

    input:
    file table from final_table

    output:
    stdout final_results

    """
    cat $table
    """
}

final_results.view{ it }
// ---------------------------------------------------------------------------------------------- //
