![Scrinvex Logo](scrinvex.png)
---
**Single Cell RNA Intron-Exon Counting**

`scrinvex` counts intronic, exonic, and junction-spanning reads for each unique barcode encountered in the input bam.
Each mapped read is checked against the input gtf to determine if the read lies entirely on introns, exons, or crosses at least one intron/exon junction.
Reads with the same UMI are only checked against any given gene once. Subsequent reads with the same UMI will not be checked against any gene that the first read intersected.

# Usage

`scrinvex {gtf} {bam} {output directory}`

# Install scR-Invex

* `git clone --recursive git@github.com:broadinstitute/scrinvex.git`
    * If you do not use `--recursive` you will be missing dependencies
    * Fix with `git submodule update --init --recursive`
* `make`

### Docker Image

scR-Invex is available via `gcr.io/broad-cga-aarong-gtex/scrinvex`

## Output Format

`scrinvex` produces 3 files in the output directory:
* `*.introns.tsv`
* `*.junctions.tsv`
* `*.exons.tsv`

Each of these three files contains 3 columns:
* `gene_id`
* `barcode`
* Counts of intron, junction, or exon reads for that gene id - barcode combination, respectively
