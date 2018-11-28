%.sorted.gff:%.gff
	gt gff3 -sort $<  > $@

%.bam: %.fa %.fa
	minimap2 -a -xasm10 $^ | samtools view -F 4 -b -o $@

%.genes.bam: %.genes.fa %.fa
	minimap2 -a -xasm10 $^ | samtools view -F 4 -b -o $@

%.reads.bam: %.fa %_R1.fq.gz %_R2.fq.gz
	bwa index $^; bwa mem $^ | samtools view -F 4 -b -o $@

%.sorted.bam: %.bam
	samtools sort $< > $@
	samtools index $@

%.: %.fa %.sorted.bam %.genes.sorted.bam %.reads.sorted.bam %.gff
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

# Extract features from GFF

%.sorted.gff: %.gff
	gt gff3 -sortlines yes -retainids yes -tidy yes -fixregionboundaries yes -addids $< > $@

%.genes.fa: %.fa %.gff
	gt extractfeat -type gene -coords -matchdescstart -retainids -seqid -seqfile $< $*-EMBOSS.gff > $@

# Table2asn
%.gbf %.sqn: %.fa %.gff
	linux64.table2asn_GFF -V bv -locus-tag-prefix Se404-851cp -i $< -f $*.gff -Z $*.dr -o $*.sqn -t $*.sbt -X C


