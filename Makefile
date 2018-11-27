%.sorted.gff:%.gff
	gt gff3 -sort $<  > $@

%.bam: %.fa %.fa
	minimap2 -a -xasm10 $^ | samtools view -F 4 -b -o $@

%.genes.bam: %.genes.fa %.fa
	minimap2 -a -xasm10 $^ | samtools view -F 4 -b -o $@

%.reads.bam: %.fa %_R1.fq.gz %_R2.fq.gz
	bwa mem $^ | samtools view -F 4 -b -o $@

%.sorted.bam: %.bam
	samtools sort $< > $@
	samtools index $@

%.: %.fa %.bam %.genes.bam %.reads.bam %.gff
	igv -g $^

%.tsv: %.genes.fa %.genes.fa
	makeblastdb -dbtype nucl -in $< -out $*
	blastn -query $*.genes.fa -out $@ -db $* -outfmt ‘6 std qcovs’ -perc_identity 100

# Using MUMmer
%.delta: %.fa
 	nucmer --maxmatch --nosimplify --prefix=$*_$* $< $<
%.coords: %.delta
 	show-coords -r $< > $@
%.repeats: %.fa
	repeat-match -n $< > $@

# Using minidot
%.pdf: %.fa
	minidot -o $@ $<

# Using minimap2
%.paf: %.fa
	minimap2 -DP $< $< > $@

%.allow-mismatches.paf: %.fa
	minimap2 -c -m200 $< $< > $@
