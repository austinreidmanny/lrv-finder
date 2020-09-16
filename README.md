# lrv-finder
A pipeline to discover Leishmaniavirus sequences from host (meta)transcriptomes

## Objective:

This pipeline will take in an accession number for a transcriptome in the NCBI Sequence Read Archive (SRA), and it will output a summary of the number of reads mapping to LRV1 or LRV2. In this way, we are looking to expand the strains of LRV1 and LRV2 and perhaps extend the range of Leishmania species that can harbor these viruses.

The pipeline is constructed as follows:
  - Download RNA-seq reads from the NBCI SRA using `fasterq-dump` from the NCBI `sratoolkit`
  - Use `BWA mem` to map any matching reads to represenatative genomes of LRV1 or LRV2.
  - Use `samtools idxstats` to translate the mapped BAM into a table of mapped counts
  - Transform the mapped .bam file to individual FASTA files per every reference sequence that has mapped reads
  - Use UNIX core utilities (cut, awk, tr ...) to format a summary of mapped reads to the user
  - Generate a final user-friendly summary to the user

------------------------------------------------------------------------------------------------

## Software:

Download the latest version of lrv-finder @ https://github.com/austinreidmanny/lrv-finder

Dependencies & requirements are:
  * Operating system: Linux/macOS
  * Nextflow (https://www.nextflow.io)
  * Conda (https://docs.conda.io/en/latest/miniconda.html)

------------------------------------------------------------------------------------------------

## Usage:

```
nextflow run lrv-finder.nf --sra "SRR123456"

optional parameters:
--output "output/dir/" [default="current-directory/output/"]
--threads INTEGER (number of CPUs) [default=4]
--memory INTEGER (amount of RAM) [in GB, default=8]
--tempdir "tmp/dir" [default="/tmp"]
--referenceseqs "reference_sequences.fasta" [default = "resources/lrv-genomes.fasta"]
```

------------------------------------------------------------------------------------------------

## Contact:

Austin R. Manny
Nibert Lab @ Harvard Medical School
austinmanny@g.harvard.edu
github.com/austinreidmanny

------------------------------------------------------------------------------------------------
