%.sorted.gff:%.gff
	gt gff3 -sort $<  > $@

