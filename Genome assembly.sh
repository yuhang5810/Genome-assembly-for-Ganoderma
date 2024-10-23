####For genome survey####
jellyfish count -C -m 31 -s 100M -t 10 <(zcat xianzhilouS1_1.fq.gz) <(zcat xianzhilouS1_2.fq.gz) -o test31S1.jf 
jellyfish histo -t 10 reads.jf > reads.histo
Rscript /pub4/yuhang/software/genomescope-master/genomescope.R reads.histo 21 150 output

####Genome assembly####
####PacBio HiFi bam2fasta####
samtools view m64064_220805_025647.subreads.bam | awk '{print ">"$1"\n"$10}' > subread.fa.gz

####hifiasm####
hifiasm --h1 hicdataR1.fq.gz --h2 hicdataR2.fq.gz -o m64064.asm -t 40 subread.fa.gz 2> test1.log

####get contigs####
gfatools gfa2fa S1.hic.p_ctg.gfa > S1.fasta
awk '{if($1=="S") print $3}' S1/S1.hic.p_ctg.gfa | awk -F "" '{print NF}' | sort -n

####ALL_HiC (option)#####
ALLHiC_pip.sh -r s2.fa -1 XZLS1_1.fq.gz -2 XZLS1_2.fq.gz -k 13 -e MboI

####purge_haplotigs####
minimap2 -ax map-pb S2.fasta S2.hifi.fa | samtools view -hF 256| samtools sort -@ 10 -m 1G -o s2.aligned.bam -T tmp.ali
purge_haplotigs readhist -b s2.aligned.bam -g S1.fasta -t 20
purge_haplotigs contigcov -i s2.aligned.bam.gencov -l 10 -m 130 -h 200
purge_haplotigs purge  -g genome.fasta  -c coverage_stats.csv  -b aligned.bam  -t 4-a 60

####ragtag####
ragtag.py scaffold genomic.fna S1.fasta -o S1/

####Pilon####
pilon --genome assembly.fasta --fix all --changes --frags assembly_illumina.sorted.bam --output pilon --outdir pilon_result  --vcf

####busco#####
busco -i S3.fasta -c 20 -m geno -l fungi_odb10 --out S3 --offline
busco -i pep.fasta -c 20 -m prot -l fungi_odb10 --out S3 --offline

####Genome annotation####
####Augusut####
perl Augustus/scripts/gff2gbSmallDNA.pl G57_HIC.gene.gff G57_HIC.fasta 1000 g57.genes.raw.db
etraining --species=generic --stopCodonExcludedFromCDS=false g57.genes.raw.db 2> train.err
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' >badgenes.lst
perl Augustus/scripts/filterGenes.pl badgenes.lst g57.genes.raw.db > genes.gb
grep '/gene' genes.gb |sort |uniq  |sed 's/\/gene=//g' |sed 's/\"//g' |awk '{print $1}' >geneSet.lst
perl getseq.pl geneSet.lst G57.pep > geneSet.lst.fa
makeblastdb -in geneSet.lst.fa -dbtype prot -parse_seqids -out geneSet.lst.fa
blastp -db geneSet.lst.fa -query geneSet.lst.fa -out geneSet.lst.fa.blastp -evalue 1e-5 -outfmt 6 -num_threads 8
cat <(awk '{if($3<70) print $1"\n"$2}' geneSet.lst.fa.blastp | sort -u) <( awk '{if($3>=70) print $1}' geneSet.lst.fa.blastp) | sort -u > gene_filter.id1
for  i in `cat gene_filter.id1`; do grep $i G57_HIC.gene.gff >> gene_filter.gff3; done
nohup perl Augustus/scripts/gff2gbSmallDNA.pl gene_filter.gff3 HIC.fasta 1000 genes.gb.filter &
perl Augustus/scripts/randomSplit.pl genes.gb.filter 100
perl Augustus/scripts/new_species.pl --AUGUSTUS_CONFIG_PATH=config --species=Ganoderma
etraining --species=spinach genes.gb.filter.train
augustus --species=spinach genes.gb.filter.test |tee firsttest.out
augustus --species=arabidopsis genes.gb.filter.test |tee firsttest_ara.out
nohup augustus --species=Ganoderma S1.fa > S1.test.gff 2>s1.out.file &

####Geta####
geta.pl --genome s1.fa \
-1 1_R1.fq.gz,2_R1.fq.gz,3_R1.fq.gz,4_R1.fq.gz\
-2 1_R2.fq.gz,2_R2.fq.gz,3_R2.fq.gz,4_R2.fq.gz\
--protein pep.fasta \
--augustus_species Ganoderma \
--RM_species Ganoderma \
--out_prefix s1 \
--cpu 20 \
--gene_prefix S1gene

####ncRNA annotation####
tRNAscan-SE -o s6tRNA.txt -f s6rRNA.ss -m s6tRNA.stast S6.fa
barrnap --threads 8 --outseq S2.rRNA.fa S2.fa > S2.rRNA.gff3
perl infernal-tblout2gff.pl --cmscan --fmt2 s2genome.tblout > s2.ncRNA.gff3
