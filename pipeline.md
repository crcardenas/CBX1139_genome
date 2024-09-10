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

## blobtoolkit, step 1
https://blobtoolkit.genomehubs.org/

references: 10.3389/fgene.2013.00237 & 10.12688/f1000research.12232.1

Blobtoolkit2 allows us to filter, and visualize the quality of our genome and potential conamination within it. We first need to create a directory for the toolkit, then we need to run additional analyses so that blobtools has data to interpret.

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

First we need to create a new directory and appropriate databases for blobtools. In general each particular database will need to be downloaded and formatted. However, this has already ben done and can be skipped. But here is an example of what it might look like. Busco, BLAST, and Diamond will each have their own database format.

Typically you would create some additional databases; but they can be found here on pyrgus due to their sizes: `/home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/data`

**See how to install blast & diamond in those sections!**

Next you will need to create a blobtools directory using blobtools:

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

BLAST is a **B**asic **L**ocal **A**lignment **S**earch **T**ool that we will use on the command line to find potential contaminants
it is a relatively efficient search tool but the database is large, so you want to provide enough resources and time to run this script. I recommend starting this first as it could take a couple days depending on the resources available on the cluster

make environment

`$ conda create -n blast -c bioconda blast`


we have already downloaded a directory for the genbank database
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

Diamond behaves a lot like blast but is generally faster, and we will use a uniprot database rather than nucleotide data base from NCBI.

Create a conda environment

`$ conda create -n diamond -c bioconda -c conda-forge diamond`

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

I concatenated the raw minion fastq files, I have already done this and you can find it here: `/home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/minimap2/raw_minion.fastq`

I would create a symlink, so you dont have to create many many big files

```
ln -s /home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/minimap2/raw_minion.fastq ./raw_minion.fastq
```

Just like before, change directories, create a bash script, and run it using nohup

```
source /local/anaconda3/bin/activate
conda activate minimap2

ASSEMLBY="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta"

minimap2 -ax map-ont \
	-t 4 ${ASSEMBLY} raw_minion.fastq | \
	samtools sort -@4 \
	-O BAM -o assembly.reads.bam
```


## busco

https://busco.ezlab.org/

references 10.1002/cpz1.323 & https://doi.org/10.1093/molbev/msab199

Busco searches for orthologous gene sequences from a standard database. This is pretty easy. Jeremy has the busco databases already setup and ready to run; so it can be run this way.

Same pattern as before, run your busco analysis and wait.

```
#! /bin/bash
source activate /home/jeremy/local/envbusco5/

GENOME="/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta"

/home/jeremy/local/envbusco5/bin/busco -i ${GENOME} \
    -o busco_${GENOME} \
    -m genome \
    -l insecta_odb10
```

## mitofinder

We need to recover mitochondrial sequences that are present in our data

```

```


## blobtoolkit, step 2

Finally, we use blob tools to evaluate our genome.

all of these outputs we've generated will be processed by blobtools. In your assembly directory run the following commands individually or as a shell script

For your blast and diamond analyses, you can, use more than one search output:
```
blobtools add \
    --hits blast/nt.out \
    --hits diamond/diamond.out \
    --taxrule bestsumorder \
    --taxdump ~/taxdump \
    CBX1139
```
or 
```
blobtools add \
    --hits blast/nt.out \
    --hits blast/nt_prok.out \
    --hits blast/nt_viruses.out \        
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

add busco scores the busco pipeline is a little clunky in its output, but this is fine, just be careful with your path here
```
blobtools add \
    --busco /home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/busco/busco_/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta/run_insecta_odb10/full_table.tsv \
    CBX1139
```

Once you have added all these tools, you can use the command line functions or use an interactive viewer. The interactive viewer will be easier to star with.


## blobtoolkit, filtering step

Using blobtools to ID contigs to remove based on our filtering critera we first run the following commmand

```
blobtools filter \
    --query-string "assembly.reads_cov--Min=20.0&length--Min=1000&bestsumorder_phylum--Keys=0%2C5%2C6%2C10%2C9%2C2%2C7%2C8%2C3%2C4&gc--Min=0.300" \
    --fasta /home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta \
    --output CBX1139_filtered \
    CBX1139
```

We can then examine this filtered output, generate blob plots, and find loci to keep based on these statistics. This will be performed later.


in the directory containing the assembly files *and* the data; use:

```
blobtools view --remote `pwd`
```

In the terminal on your __local machine__, run what the previous command tells you to do; and access the web address in your browser.

If this doesnt work, the three basic blobtoolkit plots can be created using the following command.

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
