#RNA-seq from Pseudonaja textilis milked and unmilked venom glands, including small RNA libraries. 
#TRANSCRIPTOME ASSEMBLY - P. textilis
#Reads were evaluated with Fastqc, adaptors were identified, but quality looks fine. Running Trimomatic to remove adaptors and clean up quality:
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 5 -phred33 -trimlog milked.VG.log venomous_gland_milked_RNA_ACAGTG_L002_R1_1.fastq venomous_gland_milked_RNA_ACAGTG_L002_R2_1.fastq -baseout milked.VG.filteredreads ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:30 MINLEN:30
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 5 -phred33 -trimlog unmilked.VG.log venomous_gland_unmilked_RNA_GTGAAA_L002_R1_1.fastq venomous_gland_unmilked_RNA_GTGAAA_L002_R2_1.fastq -baseout unmilked.VG.filteredreads ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:30 MINLEN:30
#ABySS assembly:
abyss-pe np=10 k=30 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_30
abyss-pe np=10 k=34 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_34
abyss-pe np=10 k=38 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_38
abyss-pe np=10 k=42 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_42
abyss-pe np=10 k=46 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_46
abyss-pe np=10 k=50 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_50
abyss-pe np=10 k=54 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_54
abyss-pe np=10 k=58 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_58
abyss-pe np=10 k=62 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_62
abyss-pe np=10 k=66 q=10 e=2 E=2 c=5 n=10 in='milked.VG.filteredreads_1P.fastq milked.VG.filteredreads_2P.fastq' name=milked.VG_66
abyss-pe np=10 k=30 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_30
abyss-pe np=10 k=34 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_34
abyss-pe np=10 k=38 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_38
abyss-pe np=10 k=42 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_42
abyss-pe np=10 k=46 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_46
abyss-pe np=10 k=50 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_50
abyss-pe np=10 k=54 q=10 e=2 E=2 c=5 n=10 in='uunmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_54
abyss-pe np=10 k=58 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_58
abyss-pe np=10 k=62 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_62
abyss-pe np=10 k=66 q=10 e=2 E=2 c=5 n=10 in='unmilked.VG.filteredreads_1P.fastq unmilked.VG.filteredreads_2P.fastq' name=unmilked.VG_66
#Merge ABySS assemblies
transabyss-merge milked.VG_30.fasta milked.VG_34.fasta milked.VG_38.fasta milked.VG_42.fasta milked.VG_46.fasta milked.VG_50.fasta milked.VG_54.fasta milked.VG_58.fasta milked.VG_62.fasta milked.VG_66.fasta --mink 30 --maxk 66
transabyss-merge unmilked.VG_30.fasta unmilked.VG_34.fasta unmilked.VG_38.fasta unmilked.VG_42.fasta unmilked.VG_46.fasta unmilked.VG_50.fasta unmilked.VG_54.fasta unmilked.VG_58.fasta unmilked.VG_62.fasta unmilked.VG_66.fasta --mink 30 --maxk 66
#Because the genome is available, a genome-guided Trinity assembly was also performed - this assembly was guided by the genome by first mapping reads, then generating contigs.
#Reads were mapped with Bowtie2 (default parameters) to the Pseudonaja textilis genome (EBS10Xv2-PRI, GCA_900518735.1). Mapping percentage was 60% for each library.
bowtie2 -p 30 --fr -x Pseudonaja_textilis_genome.bowtie2.ref -1 venomous_gland_milked_RNA_ACAGTG_L002_R1_1.fastq -2 venomous_gland_milked_RNA_ACAGTG_L002_R2_1.fastq -S Pseudonaja_textilis.milked.vg.genome.aligned.sam 2> Pseudonaja_textilis.milked.vg.genome.aligned.stats
bowtie2 -p 30 --fr -x Pseudonaja_textilis_genome.bowtie2.ref -1 venomous_gland_unmilked_RNA_GTGAAA_L002_R1_1.fastq -2 venomous_gland_unmilked_RNA_GTGAAA_L002_R2_1.fastq -S Pseudonaja_textilis.unmilked.vg.genome.aligned.sam 2> Pseudonaja_textilis.unmilked.vg.genome.aligned.stats
samtools view -bS Pseudonaja_textilis.milked.vg.genome.aligned.sam | samtools sort - Pseudonaja_textilis.milked.vg.genome.aligned.TrinityGG.ready
samtools view -bS Pseudonaja_textilis.unmilked.vg.genome.aligned.sam | samtools sort - Pseudonaja_textilis.unmilked.vg.genome.aligned.TrinityGG.ready
#Trinity genome-guided assembly:
Trinity --genome_guided_bam Pseudonaja_textilis.milked.vg.genome.aligned.TrinityGG.ready.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 30 --output Pseudonaja_textilis.milked.vg.genome.aligned.Trinity.GG.Assembly
Trinity --genome_guided_bam Pseudonaja_textilis.unmilked.vg.genome.aligned.TrinityGG.ready.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 30 --output Pseudonaja_textilis.unmilked.vg.genome.aligned.Trinity.GG.Assembly
#PEAR was used to merge paired-end reads if they overlapped, this was done separately for both milked and unmilked glands, using PEAR default parameters:
pear -f milked.VG.filteredreads_1P.fastq -r milked.VG.filteredreads_2P.fastq -o Milked.VG.mergedReads.fastq
pear -f unmilked.VG.filteredreads_1P.fastq -r unmilked.VG.filteredreads_2P.fastq -o Unmilked.VG.mergedReads.fastq
#Extender assemblies were completed with the merged reads using 10,000 starting seeds; this resulted in two separate assemblies, one for each condition.
java Extender3 Milked.VG.mergedReads.fastq Milked.VG.mergedReads.fastq 10 63 100 10000 2 10 20 Milked.VG_Extender_contigs.fasta 1 0 > extender.log &
java Extender3 Unmilked.VG.mergedReads.fastq Unmilked.VG.mergedReads.fastq 10 63 100 10000 2 10 20 Unmilked.VG_Extender_contigs.fasta 1 0 > extender.log &
#All assemblies (two Extender assemblies, two Trinity-GG assemblies and the two ABySS assemblies) were concatenated into one master fasta file
cat *.fasta > P.textilis.all.assemblies.fasta
#Remove redundancy
fastanrdb P.textilis.all.assemblies.fasta > P.textilis.vg.master.assembly.nr.fasta
#Remove transcripts less than 150 bp
cd-hit-est -i P.textilis.vg.master.assembly.nr.fasta -o P.textilis.vg.master.assembly.nr.greaterthan150bp.fasta -c 1.0 -n 5 -l 150 -T 30 -M 3000
#Organize the headers of this master file to be ready for input into the pipeline
perl /opt/evigene18may07/scripts/rnaseq/trformat.pl -output P.textilis.vg.master.assembly.nr.greaterthan150bp.ready.fasta -input P.textilis.vg.master.assembly.nr.greaterthan150bp.fasta
#Generate cds and aa sequences, best transcript is chosen from each cluster by script
/opt/evigene17mar10/scripts/prot/tr2aacds.pl -mrnaseq P.textilis.vg.master.assembly.nr.greaterthan150bp.fasta -NCPU=30 1>tr2aacds.log 2>tr2aacds.err
#Concatenate both okay and okalt transcripts
#Need to also have an abundance filter cut off (contigs with TPM < 1 removed)
#Run RSEM for Abundances
rsem-prepare-reference --bowtie2 P.textilis.VG.Master.EvidentialGene.output.complete.fasta Pseudonaja.textilis.VG.Master.rsem.ref
rsem-calculate-expression -paired-end -p 15 --bowtie2 /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/Raw_reads_and_Merged/venomous_gland_milked_RNA_ACAGTG_L002_R1_1.fastq /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/Raw_reads_and_Merged/venomous_gland_milked_RNA_ACAGTG_L002_R2_2.fastq Pseudonaja.textilis.VG.Master.rsem.ref P.textilis.VG.Master.Milked.Reads.Mapped.Apr15
rsem-calculate-expression -paired-end -p 15 --bowtie2 /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/Raw_reads_and_Merged/venomous_gland_unmilked_RNA_GTGAAA_L002_R1_1.fastq /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/Raw_reads_and_Merged/venomous_gland_unmilked_RNA_GTGAAA_L002_R2_2.fastq Pseudonaja.textilis.VG.Master.rsem.ref P.textilis.VG.Master.unMilked.Reads.Mapped.Apr15
faSomeRecords P.textilis.VG.Master.EvidentialGene.output.complete.fasta list.of.all.contigs.over1TPM.txt Final_P.textilis.VG.Master.uniq.contigs.fasta
#Final contig set above includes only contigs that are coding, unique, over 150 bp, and expressed over 1 TPM
#RNA-seq from Crotalus viridis milked venom glands, including small RNA libraries. 
#TRANSCRIPTOME ASSEMBLY - C. viridis
#Read QC with Trimomatic to remove adaptors and clean up quality:
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 5 -phred33 -trimlog CV.milked.VG.log CV.milked.VG.R1.fastq CV.milked.VG.R2.fastq -baseout CV.milked.VG.filteredreads ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:30 MINLEN:30
#Trinity de novo assembly:
Trinity --seqType fq --left CV.milked.VG.filteredreads_1P.fastq --right CV.milked.VG.filteredreads_2P.fastq --CPU 20 --max_memory 100G --output Trinity.CV.milked.VG > Trinity.CV.milked.VG.txt 
#Because the genome is available, a genome-guided Trinity assembly was also performed - this assembly was guided by the genome by first mapping reads, then generating contigs.
#Reads were mapped with Bowtie2 (default parameters) to the Crotalus viridis genome (UTA_CroVir_3.0, GCA_003400415.2). Mapping percentage was approximately 60% for the RNA-seq library.
bowtie2 -p 30 --fr -x Crotalus_viridis_genome.bowtie2.ref -1 CV.milked.VG.R1.fastq -2 CV.milked.VG.R2.fastq -S Crotalus_viridis.milked.vg.genome.aligned.sam 2> Crotalus_viridis.milked.vg.genome.aligned.stats
samtools view -bS Crotalus_viridis.milked.vg.genome.aligned.sam | samtools sort - Crotalus_viridis.milked.vg.genome.aligned.TrinityGG.ready
#Trinity genome-guided assembly:
Trinity --genome_guided_bam Crotalus_viridis.milked.vg.genome.aligned.TrinityGG.ready.bam --genome_guided_max_intron 10000 --max_memory 50G --CPU 30 --output Crotalus_viridis.milked.vg.genome.aligned.Trinity.GG.Assembly
#PEAR was used to merge paired-end reads if they overlapped, this was done separately for both milked and unmilked glands, using PEAR default parameters:
pear -f CV.milked.VG.filteredreads_1P.fastq -r CV.milked.VG.filteredreads_2P.fastq -o CV.milked.VG.mergedReads.fastq
#Extender assemblies were completed with the merged reads using 10,000 starting seeds; this resulted in two separate assemblies, one for each condition.
java Extender3 CV.milked.VG.mergedReads.fastq CV.milked.VG.mergedReads.fastq 10 63 100 10000 2 10 20 CV.milked.VG_Extender_contigs.fasta 1 0 > extender.log &
#All assemblies (two Extender assemblies, two Trinity-GG assemblies and the two ABySS assemblies) were concatenated into one master fasta file
cat *.fasta > C.viridis.all.assemblies.fasta
#Remove redundancy
fastanrdb C.viridis.all.assemblies.fasta > C.viridis.all.assemblies.nr.fasta
#Remove transcripts less than 150 bp
cd-hit-est -i C.viridis.all.assemblies.nr.fasta -o C.viridis.all.assemblies.nr.greaterthan150bp.fasta -c 1.0 -n 5 -l 150 -T 30 -M 3000
#Need to also have an abundance filter cut off (contigs with TPM < 1 removed)
#Run RSEM for Abundances
rsem-prepare-reference --bowtie2 C.viridis.all.assemblies.nr.greaterthan150bp.fasta C.viridis.VG.Master.rsem.ref
rsem-calculate-expression -paired-end -p 15 --bowtie2 CV.milked.VG.R1.fastq CV.milked.VG.R2.fastq C.viridis.VG.Master.rsem.ref  C.viridis.VG.Master.4daymilked.Reads.Mapped.Apr15
faSomeRecords C.viridis.all.assemblies.nr.greaterthan150bp.fasta list.of.all.contigs.over1TPM.txt Final_C.viridis.VG.Master.uniq.contigs.fasta
#Final contig set above includes only contigs that are unique, over 150 bp, and expressed over 1 TPM
#TRANSCRIPTOME ANNOTATIONS
#Run Diamond blastx on this final transcript subset against the NCBI nr database - for both P. textilis and C. viridis
diamond blastx -d /home/user/Documents/Databases/P_textilis_prots/P.textilis.genome.protein.set.db -q Final_P.textilis.VG.Master.uniq.contigs.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o P.textilis.final.contig.set.Blastp.genome.protein.set.matches -p 10 -k 1 --evalue 1e-5 -p 20 -b 30
diamond blastx -d /home/user/Documents/Databases/nr/Diamonddb/nr -q Final_P.textilis.VG.Master.uniq.contigs.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o P.textilis.final.contig.set.Blastp.nr.protein.set.matches -k 1 --evalue 1e-5 -p 20 -b 30
diamond blastx -d /home/user/Documents/Databases/nr/Diamonddb/nr -q Final_C.viridis.VG.Master.uniq.contigs.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o C.viridis.final.contig.set.Blastp.nr.protein.set.matches -k 1 --evalue 1e-5 -p 20 -b 30
#ORF Predictor - identifying complete toxin sequences - only P. textilis
diamond blastx -d /home/user/Documents/Databases/nr/Diamonddb/nr -q Final_P.textilis.VG.Master.uniq.contigs.fasta -o P.textilis.final.contig.set.Blastp.nr.protein.set.matches.default.output -p 30 -k 1 --evalue 1e-5 -b 30
perl /opt/OrfPredictor/OrfPredictor.pl Final_P.textilis.VG.Master.uniq.contigs.fasta P.textilis.final.contig.set.Blastp.nr.protein.set.matches.default.output 1 both e-5 Final_P.textilis.VG.Master.predicted.aa.contigs.fasta
perl /opt/OrfPredictor/extractCDS.pl Final_P.textilis.VG.Master.uniq.contigs.fasta Final_P.textilis.VG.Master.predicted.aa.contigs.fasta Final_P.textilis.VG.Master.predicted.cds.contigs.fasta
#Use of Venomix and ToxCodAn to guide toxin annotations - only P. textilis
./Venomix.py  -f /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/venomix/P.textilis_test/Final_P.textilis.VG.Master.uniq.contigs.fasta -b /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/venomix/P.textilis_test/P.textilis.ToxProt.Blast.output -d /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/venomix/P.textilis_test/P.textilis.VG.complete.toxin.contigs.Milked.Reads.Mapped.May23.genes.results
#Transcript expression levels in milked and unmilked venom glands - P. textilis
cat Final_P.textilis.VG.Master.uniq.contigs.fasta rna.fasta > Final.complete.global.plus.toxins.transcripts.fasta
rsem-prepare-reference --bowtie2 Final.complete.global.plus.toxins.transcripts.fasta Final.complete.global.plus.toxins.transcripts.rsem.ref
rsem-calculate-expression -paired-end -p 15 --bowtie2 venomous_gland_milked_RNA_ACAGTG_L002_R1_1.fastq venomous_gland_milked_RNA_ACAGTG_L002_R2_1.fastq Final.complete.global.plus.toxins.transcripts.rsem.ref Final.complete.global.plus.toxins.transcripts.Milked.Reads.Mapped.Apr15
rsem-calculate-expression -paired-end -p 15 --bowtie2 venomous_gland_unmilked_RNA_GTGAAA_L002_R1_1.fastq venomous_gland_unmilked_RNA_GTGAAA_L002_R2_1.fastq Final.complete.global.plus.toxins.transcripts.rsem.ref Final.complete.global.plus.toxins.transcripts.UnMilked.Reads.Mapped.Apr15
#Transcript expression levels in milked and unmilked venom glands
#In Rattlesnakes, downloaded from NCBI
fastq-dump --split-files SRR7401989
fastq-dump --split-files SRR7402004
fastq-dump --split-files SRR7402005
fastq-dump --split-files SRR11524062
fastq-dump --split-files SRR11524063
fastq-dump --split-files SRR11524059
fastq-dump --split-files SRR11524060
fastq-dump --split-files SRR11524050
fastq-dump --split-files SRR11524051
#In Rattlesnakes, downloaded from NCBI
rsem-prepare-reference --bowtie2 CroVir_rnd1.all.maker.transcripts.final_myotoxin.fasta C.viridis.global.plus.toxins.transcripts.bowtie2.index.rsem.ref
rsem-prepare-reference --bowtie2 GCF_016545835.1_ASM1654583v1_rna.fna C.tigris.global.transcripts.bowtie2.index.rsem.ref
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR7401989_1.fastq SRR7401989_2.fastq C.viridis.global.plus.toxins.transcripts.bowtie2.index.rsem.ref C.viridis.unmilked.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR7402004_1.fastq SRR7402004_2.fastq C.viridis.global.plus.toxins.transcripts.bowtie2.index.rsem.ref C.viridis.1daymilked.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR7402005_1.fastq SRR7402005_2.fastq C.viridis.global.plus.toxins.transcripts.bowtie2.index.rsem.ref C.viridis.3daymilked.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 CV.milked.VG.R1.fastq CV.milked.VG.R2.fastq C.viridis.global.plus.toxins.transcripts.bowtie2.index.rsem.ref C.viridis.4daymilked.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524062_1.fastq SRR11524062_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.unmilked62.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524063_1.fastq SRR11524063_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.unmilked63.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524059_1.fastq SRR11524059_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.1daymilked59.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524060_1.fastq SRR11524060_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.1daymilked60.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524050_1.fastq SRR11524050_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.4daymilked50.complete.transcriptome.plus.toxins.Reads.Mapped
rsem-calculate-expression -paired-end -p 15 --bowtie2 SRR11524051_1.fastq SRR11524051_2.fastq C.tigris.global.transcripts.bowtie2.index.rsem.ref C.tigris.4daymilked51.complete.transcriptome.plus.toxins.Reads.Mapped
#Gfold - fold change between conditions, need to prep files, use the expected counts, full length of transcript and the FPKM values to construct the input files for Gfold
/Users/cassandramodahl/Desktop/Scripts_Programs/gfold.V1.1.4/gfold diff -s1 P.textilis.unmilked.final.counts.sept21.txt -s2 P.textilis.milked.final.counts.sept21.txt -o P.text.milkedvsunmilked.Sept21.gfold.diff
/Users/cassandramodahl/Desktop/Scripts_Programs/gfold.V1.1.4/gfold diff -s1 C.viridis.unmilked.final.counts.sept20.txt -s2 C.viridis.1daymilked.final.counts.sept20.txt -o C.viridis_1day.milkedvsunmilked.Sept21.gfold.diff
/Users/cassandramodahl/Desktop/Scripts_Programs/gfold.V1.1.4/gfold diff -s1 C.viridis.unmilked.final.counts.sept20.txt -s2 C.viridis.3daymilked.final.counts.sept20.txt -o C.viridis_3day.milkedvsunmilked.Sept21.gfold.diff
/Users/cassandramodahl/Desktop/Scripts_Programs/gfold.V1.1.4/gfold diff -s1 C.tigris.unmilked62.final.counts.sept20.txt,C.tigris.unmilked63.final.counts.sept20.txt -s2 C.tigris.1daymilked59.final.counts.sept20.txt,C.tigris.1daymilked60.final.counts.sept20.txt -o C.tigris_1day.milkedvsunmilked.Sept21.gfold.diff
/Users/cassandramodahl/Desktop/Scripts_Programs/gfold.V1.1.4/gfold diff -s1 C.tigris.unmilked62.final.counts.sept20.txt,C.tigris.unmilked63.final.counts.sept20.txt -s2 C.tigris.4daymilked50.final.counts.sept20.txt,C.tigris.4daymilked51.final.counts.sept20.txt -o C.tigris_4day.milkedvsunmilked.Sept21.gfold.diff
#Human orthologs for GSEA
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Ptextilis.96.upreg.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Ptextilis.up.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Ptextilis.96.down.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Ptextilis.down.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Cviridis.96.upreg.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Cviridis.up.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Cviridis.96.down.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Cviridis.down.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Ctigris.96.up.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Ctigris.up.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
diamond blastx -d H.sapiens.uniprot.db.dmnd -q Ctigris.96.down.fasta --outfmt 6 qseqid sseqid stitle pident evalue length slen -o Ctigris.down.Human.Uniprot.orthologs.matches -k 1 --evalue 1e-5 -p 5
#Small RNA libraries 
#Fastx_cliper of the Fastx toolkit to remove adapters. 
fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 18 -i unmilkedsRNA_S2_L001_R1_001.fastq -o unmilkedsRNA.fastxclipper.fastq
fastx_clipper -a TGGAATTCTCGGGTGCCAAGG -l 18 -i milkedsRNA_S1_L001_R1_001.fastq -o milkedsRNA.fastxclipper.fastq
fastx_clipper -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -l 18 -i miRNA_Cvv_CO.fastq -o miRNA_Cvv_CO.clip.fastq
#Trimmomatic to make sure that all sequences are of high quality
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 15 -phred33 -trimlog fastxclipped.quality.unmilked.sRNA.log /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/miRNAs/Raw_reads/unmilkedsRNA.fastxclipper.fastq unmilkedsRNA.fastxclipped.quality.trimmed.fastq LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:18
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 15 -phred33 -trimlog fastxclipped.quality.milked.sRNA.log /home/user/Documents/Venomous_Snake_RNAseq/Snake_Assemblies/Summer_Pseudonaja_textilis/miRNAs/Raw_reads/milkedsRNA.fastxclipper.fastq milkedsRNA.fastxclipped.quality.trimmed.fastq LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:18
java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 5 -phred33 -trimlog fastxclipped.quality.miRNA_Cvv_CO.clip.log miRNA_Cvv_CO.clip.fastq miRNA_Cvv_CO.clipped.quality.fastq LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:18
#Removal of rRNA. Downloaded rRNA from all snakes, fasta from RNAcentral (downloaded May 3, 2019), inluding added rRNA sequences annotated from the P. textilis genome. Used Bowtie2 to make reads and all reads that mapped to rRNA were removed.
bowtie2 -p 15 -x /home/user/Documents/Databases/RNAcentral/Complete.Snake.rRNA.seqs.bowtie2.index -N 0 -U milkedsRNA.fastxclipped.quality.trimmed.fastq -S milkedsRNA.reads.alignment.snake.rRNA.sam --un norRNA.MILKED.trimmed.fastq --al rRNA 2> milkedsRNA.reads.alignment.snake.rRNA.Mapping.Stats
bowtie2 -p 15 -x /home/user/Documents/Databases/RNAcentral/Complete.Snake.rRNA.seqs.bowtie2.index -N 0 -U unmilkedsRNA.fastxclipped.quality.trimmed.fastq -S unmilkedsRNA.reads.alignment.snake.rRNA.sam --un norRNA.UNMILKED.trimmed.fastq --al rRNA_unmilked 2> unmilkedsRNA.reads.alignment.snake.rRNA.Mapping.Stats
bowtie2 -p 15 -x /home/user/Documents/Databases/RNAcentral/Complete.Snake.rRNA.seqs.bowtie2.index -N 0 -U miRNA_Cvv_CO.clipped.quality.fastq -S miRNA_Cvv_CO.clipped.quality.rRNA.sam --un norRNA.miRNA_Cvv_CO.trimmed.fastq --al rRNA_miRNA_Cvv_CO 2> miRNA_Cvv_CO.alignment.snake.rRNA.Mapping.Stats
#Shortened reads to appropriate small RNA size (18-32bp)
fastx_trimmer -l 32 -m 32 -i norRNA.MILKED.trimmed.fastq -o MILKED.norRNA.trimmed.short.fastq
fastx_trimmer -l 32 -m 32 -i norRNA.UNMILKED.trimmed.fastq -o UNMILKED.norRNA.trimmed.short.fastq
fastx_trimmer -l 32 -m 32 -i norRNA.miRNA_Cvv_CO.trimmed.fastq -o miRNA_Cvv_CO.trimmed.short.fastq
#Indexing genome for miRDeep2 and prepping files with known snake miRNAs for homolog analysis
bowtie-build P.textilis.genome.fasta P.textilis.genome.bowtie1.ref
bowtie-build GCA_003400415.2_UTA_CroVir_3.0_genomic.fna C.vv.genome.bowtie1.ref
remove_white_space_in_id.pl P.textilis.genome.fasta > P.textilis.genome.fa
remove_white_space_in_id.pl GCA_003400415.2_UTA_CroVir_3.0_genomic.fna > Cvv.genome.ready.fa
remove_white_space_in_id.pl miRNA_Cvv_CO.trimmed.short.fasta > miRNA_Cvv_CO.trimmed.short.ready.fasta
remove_white_space_in_id.pl known.snake.miRNAs.fastaa > known.snake.miRNAs.ready.fasta
/opt/mirdeep2-master/bin/rna2dna.pl known.snake.miRNAs.ready.fasta > all.known.snake.miRNAs.ready.fa
#Formatting reads for mapper
fastq_to_fasta -i MILKED.norRNA.trimmed.short.fastq -o MILKED.norRNA.trimmed.short.fasta
fastq_to_fasta -i UNMILKED.norRNA.trimmed.short.fastq -o UNMILKED.norRNA.trimmed.short.fasta
fastq_to_fasta -i miRNA_Cvv_CO.trimmed.short.fastq -o miRNA_Cvv_CO.trimmed.short.fasta
#Mapping reads to genomes
mapper.pl MILKED.norRNA.trimmed.short.rename.fasta -c -j -m -p P.textilis.genome.bowtie1.ref -s MILKED.miRNA.collapsed.reads.fa -t MILKED.miRNA.collapsed.reads.genome.mapped.arf -v
mapper.pl UNMILKED.norRNA.trimmed.short.rename.fasta -c -j -m -p P.textilis.genome.bowtie1.ref -s UNMILKED.miRNA.collapsed.reads.fa -t UNMILKED.miRNA.collapsed.reads.genome.mapped.arf -v
mapper.pl miRNA_Cvv_CO.trimmed.short.ready.fasta -c -j -m -p C.vv.genome.bowtie1.ref -s Cvv.CO.miRNA.collapsed.reads.fa -t Cvv.CO.miRNA.collapsed.reads.genome.mapped.arf -v
#Run miRDeep2 to identify miRNAs and quantify
miRDeep2.pl MILKED.miRNA.collapsed.reads.fa P.textilis.genome.fa MILKED.miRNA.collapsed.reads.genome.mapped.arf none Ophiophagus_hannah_AND_rna_typemiRNA.ready.fa none 2> milked.miRNA.miRDeep2.report.log
miRDeep2.pl UNMILKED.miRNA.collapsed.reads.fa P.textilis.genome.fa UNMILKED.miRNA.collapsed.reads.genome.mapped.arf none Ophiophagus_hannah_AND_rna_typemiRNA.ready.fa none 2> unmilked.miRNA.miRDeep2.report.log
miRDeep2.pl Cvv.CO.miRNA.collapsed.reads.fa Cvv.genome.ready.fa Cvv.CO.miRNA.collapsed.reads.genome.mapped.arf none all.known.snake.miRNAs.ready.fa none 2> Cvv.CO.miRNA.miRDeep2.report.log
#Target prediction with miRanda
#Toxin prediction for all miRNAs over 100 CPM
miranda Cvv.miRNA.over100CPM.fasta CVV.genome.annotated.toxin.seqs.fasta -en -19 > Cvv.100CPM.milked.miRNA.TOXINS.minus19targets.out
miranda Milked.miRNA.100CPM.miRNAs.fasta P.textilis.milked.study.TOXINS.fasta -en -19 > P.text.100CPM.milked.miRNA.TOXINS.minus19targets.out
miranda Unmilked.miRNA.100CPM.miRNAs.fasta P.textilis.milked.study.TOXINS.fasta -en -19 > P.text.100CPM.unmilked.miRNA.TOXINS.minus19targets.out
#Global transcriptome miRNA target predictions
miranda Cvv.miRNA.over100CPM.fasta CVV.genome.annotated.toxin.seqs.fasta -en -30 > Cvv.100CPM.milked.miRNA.TOXINS.minus30targets.out
miranda Milked.miRNA.100CPM.miRNAs.fasta P.textilis.milked.study.TOXINS.fasta -en -30 > P.text.100CPM.milked.miRNA.TOXINS.minus30targets.out
miranda Unmilked.miRNA.100CPM.miRNAs.fasta P.textilis.milked.study.TOXINS.fasta -en -30 > P.text.100CPM.unmilked.miRNA.TOXINS.minus30targets.out
grep -i -B 10 'Energy:' Cvv.100CPM.milked.miRNA.10TPM.minus30targets.out > Cvv.100CPM.milked.miRNA.10TPM.minus30targets.final.out
grep -i -B 10 'Energy:' P.text.100CPM.milked.miRNA.10TPM.minus30targets.out > P.text.100CPM.milked.miRNA.10TPM.minus30targets.final.out
grep -i -B 10 'Energy:' P.text.100CPM.unmilked.miRNA.10TPM.minus30targets.out > P.text.100CPM.unmilked.miRNA.10TPM.minus30targets.final.out