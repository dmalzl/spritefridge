# spritefridge
A python toolbox for processing SPRITEseq data

## Installation
To be able to run everything correctly we need a few prerequisits installed especially bedtools. 
Furthermore, at the time of writing this some dependencies refused to compile when installing with pip (`krbalancing`).
Installing these is easiest done using conda. For convenience we provide an environment file (`env.yml`) with this package
Installation thus works like
```
conda env create -f env.yml
conda activate sprite
pip install spritefridge
```

## Usage
`spritefridge` comprises five tools to process and annotate SPRITE-seq data and results. Below are some example commands. For more details please
refer to the generated help messages `spritefridge <subcommand> -h`

### extractbc
`extractbc` aims to extract barcodes from reads according to a list of used barcodes and barcode layouts (i.e. how the barcodes are aranged in read sequence)
A typical command looks like this
```
spritefridge extractbc \
    -r1 r1.fq.gz \
    -r2 r2.fq.gz \
    -bc barcodes.tsv \
    -l1 DPM \
    -l2 'Y|SPACER|ODD|SPACER|EVEN|SPACER|ODD' \
    -m 'DPM:0,Y:0,EVEN:2,ODD:2' \
    -o out.bcextract.fq.gz \
    -p 4
```
This command will read in the barcodes and the try to find barcodes in the respective read sequence in the order given by the layouts starting from 5' end.
`-m` gives the allowed mismatches for the barcode identification. In addition to `out.bcextract.fq.gz` which contains reads with the extracted barcodes appended to their names, the tool also outputs statistics for how many reads were found with 1, 2, 3, ... barcodes. `-p` specifies the number of processes to use for extraction.

### pairs
`pairs` identifies barcode clusters from aligned reads and writes them into pairs files for each cluster size
```
spritefridge pairs \
    -b in.bam \
    -o pairs/out \
    -cl 2 \
    -ch 1000 \
    --separator '['
```
This command will read in alignments from `in.bam` (needs to be filtered for multimappers and quality) groups the reads by barcodes and then writes all possible pairs for each cluster of sizes between 2 and 1000 reads to a file named `pairs/out_<clustersize>.pairs`. This tool also outputs a dedicated bedfile containing all reads from each cluster to be used to annotated the Cooler bins later on (see `annotate`)

### combine
`combine` merges cool files generated from cluster pairs files according to the SPRITE-seq recommendation by multiplying the counts of each Cooler by 2/n,
where n is the cluster size, before merging. The cluster size is inferred from the file name which needs to be of the pattern `<name>_<clustersize>.cool`
```
spritefridge combine \
    -i coolers/* \
    -o merged.cool \
    --floatcounts
```
`--floatcounts` ensures that merged counts are stored as float and not be casted to int

### annotate
`annotate` takes in a bedfile (see `pairs`) and annotated each bin with the overlapping reads of each cluster.
```
spritefridge annotate \
    -i merged.mcool \
    -b clusters.bed
```
`merged.mcool` is a zoomified version of the `merged.cool` file

### balance
`balance` is used to balance the contact matrices of the resulting mcool file using iterative correction and Knight-Ruiz matrix balancing
genomewide and per chromosome
```
spritefridge balance \
    -m testdata/sprite.new.mcool \
    -p 2 \
    --overwrite
```
`-p` specifies the number of processes to use for iterative correction and `--overwrite` will overwrite any existing weights with the same name in the Cooler
