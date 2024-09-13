To get genomes, use:
`datasets download genome accession --inputfile ../genome_accessions.txt --assembly-level chromosome --include genome`

simlink the fna files to here
`for i in ncbi_dataset/data/*/*.fna; do FNA=$(echo ${i} | cut -d "/" -f 5); ln -s ${PWD}/${i} ${PWD}/${FNA}; done`

running busco like:
```
#! /bin/bash

#shell script to run busco
# run like: bash busco.sh genomic.fna
# or for a directory with FNA files: nohup bash -c 'for i in 1_genomes/*.fna; do bash busco.sh ${i}; done' > busco.log &
source activate /home/jeremy/local/envbusco5/

GENOME=$1

/home/jeremy/local/envbusco5/bin/busco -i ${GENOME} \
	-o ./2_busco/busco_${GENOME} \
	-m genome \
	-l insecta_odb10
```
once job is complete, need to append sample ID fasta header and concatenate shared busco genes

make symlinks

```
for i in ${PWD}/2_busco/busco_1_genomes/*; do DIRECT=$(echo ${i} | cut -d "/" -f  9 | cut -d "." -f 1); echo ln -s ${i}/run_insecta_odb10/busco_sequences/single_copy_busco_sequences ${PWD}/3_matrices/${DIRECT}; done
```

To get the CBX genome to play nice with this pipeline, append GCA to a symlink (GCA_CBX1139). *a hacky solution*

```
ln -s /home/cody/CALOSOMA_Genomes/NANOPORE/blobtools/CBX1139/busco/busco_/home/cody/CALOSOMA_Genomes/NANOPORE/data/CBX1139_flye_test/assembly.fasta/run_insecta_odb10/busco_sequences/single_copy_busco_sequences ${PWD}/GCA_CBX1139
```

#### List of busco's to ignore from CBX1139; they map to contaminated reads. So they can just be removed from the phylogenetic analysis
54419at50557 70122at50557 144087at50557 19754at50557 84610at50557 129072at50557 143030at50557 30467at50557 47269at50557 64673at50557 80996at50557 123703at50557 149245at50557 12466at50557 70279at50557 70721at50557 108777at50557 129103at50557 144668at50557

So, we use the following command in the  CBX1139 directory
```
for i in 54419at50557 70122at50557 144087at50557 19754at50557 84610at50557 129072at50557 143030at50557 30467at50557 47269at50557 64673at50557 80996at50557 123703at50557 149245at50557 12466at50557 70279at50557 70721at50557 108777at50557 129103at50557 144668at50557; do rm ${i}.fna; done
```

concatenate and clean up headers
```
# first creat temporary directory
for i in GCA*/; do
NAME=$(echo ${i} | cut -d "/" -f 1);
mkdir tmp_${NAME}; cp ${NAME}/*.fna tmp_${NAME}/;
cd tmp_${NAME};
	for f in *.fna; do 
	FASTA=$(echo ${f} | cut -d "." -f 1);
	sed -i "s/^>.*/>${f%.*}/" ${f};
	awk '!/^>/ {printf "%s", $0; n="\n"}; /^>/ {print n $0; n=""}; END {printf "%s", n}' ${f} > ${FASTA}.tmp;
	done;
cd ..;
cat tmp_${NAME}/*.tmp > ${NAME}.tmp;
rm -r tmp_${NAME};
done
```

make a new directory, `loci` and run the following command from within that directory, to create multi fastas for alignment

```
awk '/^>/{locus_name=$1;gsub(">","",locus_name); next}{sequence=$0;sample_name=FILENAME; gsub(".tmp", "", sample_name); gsub("../", "",sample_name); output_file=locus_name".fasta";print ">"sample_name >> output_file;print sequence >> output_file}' ../*.tmp
```

Align loci now

```
#!/bin/bash
# activate conda environment
source /local/anaconda3/bin/activate
conda activate mafft
for i in loci/*.fasta; 
    do NAME=$(echo ${i} | cut -d "/" -f 2 | cut -d "." -f 1); 
    mafft --thread 3 --auto loci/${NAME}.fasta > \
    alignments/${NAME}.aligned.fasta; 
done
```

use trimal to clean up
```
#!/bin/bash
# activate conda environment
source /local/anaconda3/bin/activate
conda activate trimal
for i in alignments/*.aligned.fasta 
    do NAME=$(echo ${i} | cut -d "/" -f 2 | cut -d "." -f 1); 
    trimal -automated1 -in alignments/${NAME}.aligned.fasta -fasta \
    -out trimmed/${NAME}.trimmed.fasta; 
done
```

OK, now use amas to concatenate the loci together

```
~/AMAS/amas/AMAS.py concat -i trimmed/*.fasta -f fasta -d dna -u nexus
```

Once matrix is generated use iqtree2 to generate partitioning and then Sh-ALRT and UFBoot trees
```
#!/bin/bash

cd /home/users/c/cardenac/Cdepresus_genome_phylo

sbatch --job-name C_depressus_genome \
--mem=56000 \
--ntasks 1 \
--cpus-per-task 14 \
--partition public-cpu \
--time 3-00:00:00 \
--wrap "
module load Anaconda3;
source activate iqtree2;
iqtree -s concatenated.nex \
-p parti.nex \
-m MF+MERGE \
-rclusterf 10 \
-T 14 \
-safe;
```

```
#!/bin/bash

cd /home/users/c/cardenac/Cdepresus_genome_phylo

sbatch --job-name C_depressus_genome_UFB \
--mem=40000 \
--ntasks 1 \
--cpus-per-task 9 \
--partition public-cpu \
--time 2-00:00:00 \
--wrap "
module load Anaconda3;
source activate iqtree2;
iqtree -s concatenated.nex -p parti.nex.best_model.nex \
-o GCA_963971575 -allnni -bnni -bb 1000 -alrt 1000 \
-T 9 -safe ;
"
```
