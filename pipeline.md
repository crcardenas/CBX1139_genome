# Carabus (Platycarabus) depresus (CBX1139) genome pipeline 

Ethanol preserved specimen collected at ITALY: Valle d'Aosta, Saint-Rhémy-en-Bosses, Col du Grand Saint-Bernard in 2019.
Performed a simple high molecular weight DNA Extraction using DNA blood and tissue kit reagents.

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

Run using flye's automatic parameter selection
	
```
flye --nano-corr [PATH] --out-dir [PATH] --threads 4
```

flye found duplicate sequences, likely due to having to restart dorado.

Removing duplicates using seqkit:
```
cat CBX1139_processed/CBX1139.fastq | seqkit rmdup -n -D CBX1139.fq.duplicated.details.txt -o CBX1139_processed/CBX1139.rmdup.fastq
```

rerunning flye like so (use nohup so the job runs in the background)
```
flye --nano-corr [PATH] --out-dir [PATH] --threads 6
```

## blobtoolkit, step 1
https://blobtoolkit.genomehubs.org/

references: 10.3389/fgene.2013.00237 & 10.12688/f1000research.12232.1

```
conda create -n btk -c conda-forge python=3.9
conda activate btk
pip install "blobtoolkit[full]"
```

Create a blobtools directory using blobtools:
```
blobtools create --fasta [path/to/your/assembly.fasta] CBX1139
```

Add meta data for our sample, find the meta.json file (`./CBX1139/meta.json`) the field of interest will look like this:

```
...
"taxon": {
    "phylum":"Arthropoda",
    "class":"Insecta",
    "order":"Coleoptera",
    "suborder":"Adephaga",
    "family":"Carabidae",
    "subfamily":"Carabinae",
    "tribe":"Carabin",
    "genus":"Carabus",
    "subgenus":"Platycarabus",
    "species":"depressus",
    "name":"Carabus depressus",
    "species name":"Carabus depressus",
    "kingdom":"Metazoa"
},
...
```

## BLAST

```
#!/bin/bash
source /local/anaconda3/bin/activate
conda activate blastn

for i in nt nt_others nt_prok nt_viruses; do
        ASSEMLBY="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta";
        blastn -db ${i} -query ${ASSEMLBY} \
        -outfmt "6 qseqid staxids bitscore std" \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -evalue 1e-25 \
        -num_threads 4 \
        -out ${i}.out;
done
```

then run

```
nohup bash multi-bash.sh > multi-bash.out &!
```

## diamond
https://github.com/bbuchfink/diamond

reference doi:10.1038/nmeth.3176

Use the uniprot database rather than nucleotide data base from NCBI.

```
source /local/anaconda3/bin/activate
conda activate diamond

ASSEMLBY="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta"

diamond blastx \
        --query ${ASSEMLBY} \
        --db ../../data/uniprot/reference_proteomes.dmnd \
        --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --sensitive \
        --max-target-seqs 1 \
        --evalue 1e-25 \
        --threads 4 \
        --verbose \
        > diamond.out
```

## minimap2
https://github.com/lh3/minimap2

references doi:10.1093/bioinformatics/bty191 & doi:10.1093/bioinformatics/btab705

`$ conda create -n minimap2 -c bioconda minimap2`

I concatenated the raw minion fastq files for mapping: `/home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/minimap2/raw_minion.fastq`

```
source /local/anaconda3/bin/activate
conda activate minimap2

ASSEMLBY="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta"

minimap2 -ax map-ont \
	-t 4 ${ASSEMBLY} raw_minion.fastq | \
	samtools sort -@4 \
	-O BAM -o assembly.reads.bam
```

Then use awk to calculate average depth:
```
samtools depth assembly.reads.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

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

## search for mtDNA contigs

Recover mitochondrial sequences that are present in our data

```
# waiting on baobab

# script will show how mtDNA_contigs_keep.list (see end of MD file) was produced
```

## blobtoolkit, step 2

Use blob tools to evaluate our genome.

all of these outputs we've generated will be processed by blobtools. In your assembly directory run the following commands individually or as a shell script

For your blast and diamond analyses use more than one search output:
```
blobtools add \
    --hits blast/nt.out \
    --hits diamond/diamond.out \
    --taxrule bestsumorder \
    --taxdump ~/taxdump \
    CBX1139
```

Add coverage:

```
blobtools add \
    --cov minimap2/assembly.reads.bam \
    CBX1139
```

add busco scores:
```
blobtools add \
    --busco /home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/busco/busco_/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta/run_insecta_odb10/full_table.tsv \
    CBX1139
```
## blobtoolkit, filtering step

Using blobtools to ID contigs to remove based on our filtering critera we first run the following commmand

```
blobtools filter \
    --query-string "assembly.reads_cov--Min=20.0&length--Min=1000&bestsumorder_phylum--Keys=0%2C5%2C6%2C10%2C9%2C2%2C7%2C8%2C3%2C4&gc--Min=0.300" \
    --fasta /home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta \
    --output CBX1139_filtered \
    CBX1139
```

generate basic descriptive blobtools plots:

```
for PLOT in blobl cumulative snail; do
	blobtools view --plot \
		--format svg \
		--view ${PLOT} \
		--out ./blobtools/CBX1139/ \
		./blobtools/CBX1139/;
done
```

Last, we use seqkit to remove contigs from the assembly files using:
```
seqkit grep -f contigs_keep.list  /home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta -o nuc_DNA_CBX1139_filtered.fasta
seqkit grep -f mtDNA_contigs_keep.list /home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta -o mit_DNA_CBX1139_filtered.fasta
```
