#!/bin/bash
IFS=$'\n'
#file=$1
file="Se404-851cp.gff"
temp="temp.gff"
output="Se404-851cp.numbered.gff"
gene_count=1
mRNA_count=1
tRNA_count=1
rRNA_count=1
code="Se404-851-"
feat="gene"

if [ -e $output ]
then
	rm $output
fi

sed 's/Parent=/ID=/' $file > $temp
sed -i 's/biological_region/mRNA/' $temp
sed -i 's/[[:space:]]\+$//' $temp
sed -i 's/$/;/' $temp

for line in $(cat $temp)
do
	feature=$(echo $line | awk '{print $3}')
	case $feature in
		gene)
			name=$(echo $line | awk '/gene=/ {print $9}' | awk -F ";" '{print $2}' | awk -F "=" '{print $2}')
			echo $line | sed "s/ID=${code}\.[0-9]\+/ID=gene${gene_count}/" | sed "s/$/Name=${name}/" >> $output
			gene_count=$((gene_count + 1)) 
			feat="gene"
			;;
		mRNA)
			echo $line |  sed "s/ID=${code}\.[0-9]\+/ID=mRNA${mRNA_count}/" | sed "s/gene=\w\+-\?\w*;\?/Parent=gene$((gene_count-1));/" | sed "s/$/gene=${name}/" | sed 's/\t[0-2]\tID/\t\.\tID/' >> $output
			mRNA_count=$((mRNA_count+1))
			feat="mRNA"
			;;
		rRNA)
			echo $line | sed "s/ID=${code}\.[0-9]\+/ID=rRNA${rRNA_count}/" | sed "s/gene=\w\+-\?\w*;\?/Parent=gene$((gene_count-1));/" | sed "s/$/gene=${name}/" >> $output
			rRNA_count=$((rRNA_count+1))
			feat="rRNA"
			;;
		tRNA)
			echo $line | sed "s/ID=${code}\.[0-9]\+/ID=tRNA${tRNA_count}/" | sed "s/gene=\w\+-\?\w*;\?/Parent=gene$((gene_count-1));/" | sed "s/$/gene=${name}/" >> $output
			tRNA_count=$((tRNA_count+1))
			feat="tRNA"
			;;
		exon)
			
			case $feat in
				mRNA) 
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=mRNA$((mRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/"  >> $output
					;;
				tRNA)
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=tRNA$((tRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/"  >> $output
					;;
				rRNA)
					echo $line |  sed "s/ID=${code}\.[0-9]\+/Parent=rRNA$((rRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//'| sed "s/$/gene=${name}/"  >> $output
					;;
				*)
					;;
			esac
			;;
		CDS)
			case $feat in
				mRNA) 
					echo $line |  sed "s/ID=${code}\.[0-9]\+/Parent=mRNA$((mRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/" | sed 's/$/;/' |  sed "s/$/Name=CDS-${name}/" >> $output
					;;
				tRNA)
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=tRNA$((tRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/" | sed 's/$/;/' | sed "s/$/Name=CDS-${name}/"  >> $output
					;;
				rRNA)
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=rRNA$((rRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/" | sed 's/$/;/' | sed "s/$/Name=CDS-${name}/"  >> $output
					;;
				*)
					;;
			esac
			;;
		intron)
			case $feat in
				mRNA) 
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=mRNA$((mRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/"  >> $output
					;;
				tRNA)
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=tRNA$((tRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/"  >> $output
					;;
				rRNA)
					echo $line | sed "s/ID=${code}\.[0-9]\+/Parent=rRNA$((rRNA_count-1))/" | sed 's/gene=\w\+-\?\w*;\?//' | sed "s/$/gene=${name}/"  >> $output
					;;
				*)
					;;
			esac
			;;
		*)
			echo $line >> $output
			;;
esac
done
sed -i 's/Se404-851-/Se404-851cp/' $output
sed -i 's/featflags=type:CDS;//' $output
sed -i 's/number=[0-9];\?//' $output
sed -i 's/;$//' $output
rm $temp
