![Scrinvex Logo](scrinvex.png)
---
**Single Cell RNA Intron-Exon Counting**


`scrinvex` counts intronic, exonic, and junction-spanning reads for each unique barcode encountered in the input bam.
Each mapped read is checked against the input gtf to determine if the read lies entirely on introns, exons, or crosses at least one intron/exon junction.
Reads with the same UMI are only checked against any given gene once. Subsequent reads with the same UMI will not be checked against any gene that the first read intersected.

# Usage

`scrinvex {gtf} {bam} {output directory}`

## Output Format

`scrinvex` produces 3 files in the output directory:
* `*.introns.tsv`
* `*.junctions.tsv`
* `*.exons.tsv`

Each of these three files contains 3 columns:
* `gene_id`
* `barcode`
* Counts of intron, junction, or exon reads for that gene id - barcode combination, respectively
