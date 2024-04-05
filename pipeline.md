# Carabus (Platycarabus) depresus (CBX1139) genome pipeline 

Ethanol preserved specimen collected at ITALY: Valle d'Aosta, Saint-Rhémy-en-Bosses, Col du Grand Saint-Bernard in 2019.
Performed a simple high molecular weight DNA Extraction using DNA blood and tissue kit.

### DNA extraction QC results
Fragment analysis~ 24.3 Kb peak; ~ 21.9 Kb smear (based on 2.5 ng/uL dilution). Elution volume of 100uL ddH2O; first qubit 1.085 μg; post-cleanup qubit 6.04 ng/μL or 543.6 ng total weight. Molarity for sequencing: 40.3 fmol

## genome assembly pipeline overview
1. dorado to perform base calling
2. assembly using flye (just contigs, no scaffolds)
3. blobtoolkit
4. blast & diamond search for contamination
5. minimap2; get sequencing depth information from raw sequences
6. busco information
7. filtering & concatenation

# Running scripts

on pyrgus... nohup bash myscript.sh > myscript.out &!
or 
nohup bash myscript.sh &! 

the second will produce a nohup.out file, but will be overwritten every time.

## Dorado
https://github.com/nanoporetech/dorado

Using the raw reads from the minion we reran the base calling using a GPU and dorado

- `dorado_0.3.1_8218837_CBX1139_basecalling-res.bam` is the final dorado output. 
- `dorado_0.3.1_*.bam` indicates the dorado version (this is what was installed on baobab)
- `*_8218837_*.bam` is the slurm job
- `*_CBX1139_basecalling*.bam` is the sample and dorado program 	ran
- `*-res.bam` indicates the basecalling was rerun

First sorted bam file using samtools, then I used bedtools to convert it to fastq. Final output is CBX1139.fastq. Retrieved basic statistics using seqkit stat -Ta 


## Flye
https://github.com/fenderglass/Flye

reference doi:10.1038/s41587-019-0072-8

Fly is a genomic assembler for metagenomic or single genome assemblies. 

Create a conda environment for flye like so:
```
$ conda create -n flye -c bioconda flye
$ conda activate flye
```

to run, it should be fairly straight forward, using flye's 	automatic parameter selection
	
```
flye --nano-corr [PATH] --out-dir [PATH] --threads 4
```

flye found duplicate sequences, likely due to having to restart dorado.

Removing duplicates using:
```
cat CBX1139_processed/CBX1139.fastq | seqkit rmdup -n -D CBX1139.fq.duplicated.details.txt -o CBX1139_processed/CBX1139.rmdup.fastq
```
(you will need to create an environment that has seqkit like so: `$ conda create -n seqkit -c bioconda seqkit` and activating that environment)

rerunning flye like so (use nohup so the job runs in the background)
```
flye --nano-corr [PATH] --out-dir [PATH] --threads 6
```

## blobtoolkit
https://blobtoolkit.genomehubs.org/

references: 10.3389/fgene.2013.00237 & 10.12688/f1000research.12232.1

blobtools is not entirely intuitive, I will work on filling out this section as I wrap my head around my understanding of it


install blobtoolkit

```
conda create -n btk -c conda-forge python=3.9
conda activate btk
pip install "blobtoolkit[full]"
```
(pip is a different package manager, don't worry too much about this)

check that blobtoolkit is installed by running `$ blobtools -h`; you should see something along the lines of:

```
blobtools

BlobTools2 - assembly exploration, QC and filtering.

usage: blobtools [<command>] [<args>...] [-h|--help] [--version]
...
```

First we need to create a new directory and appropriate databases for blobtools. In an appropriate directory run the following
```
mkdir -p taxdump;
cd taxdump;
curl -L ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | tar xzf -;
cd ..;
```
**See how to install blast & diamond in those sections!**

Typically you would create some additional databases; but they can be found here on pyrgus due to their sizes: `/data/raw/Calosoma/NANOPORE/blobtools/data/blastdb` and `/data/raw/Calosoma/NANOPORE/blobtools/data/uniprot`

Next you will need to create a blobtools directory using blobtools:

```
blobtools create --fasta [path/to/your/assembly.fasta] CBX1139
```

you will need to create directories for each analysis you want to run! It helps organize the data and keep things sane (trust me)

Your diectory should look something like this:

```
./CBX1139/
├── blasts
├── busco
├── diamond
├── gc.json
├── identifiers.json
├── length.json
├── meta.json
├── minimap2
├── multi-blash.out
├── multi-blash.sh
└── ncount.json
```
## BLAST

BLAST is a **B**asic **L**ocal **A**lignment **S**earch **T**ool that we will use on the command line to find potential contaminants
it is a relatively efficient search tool but the database is large, so you want to provide enough resources and time to run this script. I recommend starting this first as it could take a couple days depending on the resources you have available

make environment

`$ conda create -n blast -c bioconda blast`


we have already downloaded a directory for the genbank database (see blobtools section)

```
blastn -db /data/raw/Calosoma/NANOPORE/blobtools/data/blastdb/nt \
    -query [/path/to/your/assembly.fasta] \
    -outfmt "6 qseqid staxids bitscore std" \
    -max_target_seqs 10 \
    -max_hsps 1 \
    -evalue 1e-25 \
    -num_threads 4 \
    -out blast.out
```

## diamond
https://github.com/bbuchfink/diamond

reference doi:10.1038/nmeth.3176

`$ conda create -n diamond -c bioconda -c conda-forge diamond`

```
diamond blastx --query [/path/to/your/assembly.fasta] \
    --db /data/raw/Calosoma/NANOPORE/blobtools/data/uniprot/uniprot.db.with.taxids \
    --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --sensitive \
    --max-target-seqs 1 \
    --evalue 1e-25 \
    --threads 4 \
    > diamond.out
```

## minimap2
https://github.com/lh3/minimap2

references doi:10.1093/bioinformatics/bty191 & doi:10.1093/bioinformatics/btab705

`$ conda create -n minimap2 -c bioconda minimap2`



## busco

https://busco.ezlab.org/

references 10.1002/cpz1.323 & https://doi.org/10.1093/molbev/msab199

```
#! /bin/bash
source activate /home/jeremy/local/envbusco5/

GENOME="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta"

/home/jeremy/local/envbusco5/bin/busco -i ${GENOME} \
    -o busco_${GENOME} \
    -m genome \
    -l insecta_odb10
```

## blobtools

