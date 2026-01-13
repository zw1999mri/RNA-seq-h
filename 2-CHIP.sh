
jun=MG1655 
mkdir 1-QC
mkdir 2-filter-data
mkdir 3-bam
mkdir 4-count


ln -s /home/zhangwei/zw/code/ref/${jun}/${jun}.gtf
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa .
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa.amb .
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa.ann .
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa.bwt .
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa.pac .
cp /home/zhangwei/zw/code/ref/${jun}/${jun}.fa.sa .




ls 0-rawdata/*R1.fastq.gz >1
ls 0-rawdata/*R2.fastq.gz >2
paste 1 2
paste 1 2 >config
cat config 

mkdir 2-filter-data
dir=./2-filter-data
cat config |while read id
do
    arr=($id)
    fq1=${arr[0]}
    fq2=${arr[1]}
fastp \
	-i $fq1 \
	-I $fq2 \
	-o ${fq1%R1.fastq.gz}\val_1.fq.gz \
	-O ${fq2%R2.fastq.gz}\val_2.fq.gz \
	--cut_front 20 \
	--cut_by_quality3 20 \
	--qualified_quality_phred 30 \
	--unqualified_percent_limit 20 \
	--n_base_limit 5 \
	--length_required 50 \
	-h ${fq1%_R1.fastq.gz}.html \
	--thread 8
done



mv 0-rawdata/*.html 1-QC/
mv 0-rawdata/*val_1.fq.gz 2-filter-data/
mv 0-rawdata/*val_2.fq.gz 2-filter-data/

ls 1-QC|while read id; do echo $id >>QC.txt; grep 'reads passed filters' 1-QC/${id} | awk -F'[()]' '{for(i=2;i<=NF;i+=2) print $i}'>>QC.txt; done


ls 2-filter-data/*val_1.fq.gz >3
ls 2-filter-data/*val_2.fq.gz >4
paste 3 4
paste 3 4  >filter-condig
cat filter-condig 
cat filter-condig |while read idd
do
    arr=($idd)
    fq1=${arr[0]}
    fq2=${arr[1]}
bwa mem -t 30 ${jun}.fa $fq1 $fq2>$fq1.sam
done

ls 2-filter-data/*.sam |while read id;do (samtools sort -O bam -@ 5 -o 3-bam/$(basename ${id} "_val_1.fq.gz.sam").bam ${id});done




#ls 3-sam |while read idd
#do 
#	samtools sort 3-sam/$idd -o 3-bam/$idd.sorted.bam
#done

t1=a
t2=b
t3=c
t4=d
t5=e
t6=f
c1=g
c2=h
c3=i
c4=j
c5=g
c6=h

macs2 callpeak -t 3-bam/${t1}.bam -c 3-bam/${c1}.bam -f BAM -g 46416520 -n WA -B -q 0.001 --outdir macs-result
macs2 callpeak -t 3-bam/${t2}.bam -c 3-bam/${c2}.bam -f BAM -g 46416520 -n WB -B -q 0.001 --outdir macs-result
macs2 callpeak -t 3-bam/${t3}.bam -c 3-bam/${c3}.bam -f BAM -g 46416520 -n QA -B -q 0.001 --outdir macs-result
macs2 callpeak -t 3-bam/${t4}.bam -c 3-bam/${c4}.bam -f BAM -g 46416520 -n QB -B -q 0.001 --outdir macs-result
macs2 callpeak -t 3-bam/${t5}.bam -c 3-bam/${c5}.bam -f BAM -g 46416520 -n WC -B -q 0.001 --outdir macs-result
macs2 callpeak -t 3-bam/${t6}.bam -c 3-bam/${c6}.bam -f BAM -g 46416520 -n WD -B -q 0.001 --outdir macs-result


ls macs-result/*peaks.xls|while read id
do
echo $id
idd=$(basename $id "_peaks.xls")
awk '{print $1"\t"$2"\t"$3"\t"$9"\t"$10}'  $id|grep '_peak_' > $id.bed
bedtools getfasta -fi MG1655.fa -bed $id.bed -fo $id.fa
done



ls macs-result/*peaks.xls|while read id
do
echo $id
idd=$(basename $id "_peaks.xls")
bedtools getfasta -fi MG1655.fa -bed $id.bed -fo $id.fa
done


mkdir 3-meme

ls macs-result/*peaks.xls|while read id
do
echo $id
idd=$(basename $id "_peaks.xls")
meme $id.fa -dna -oc ${id}meme -nostatus -time 1440000 -mod zoops  -minw 6 -maxw 20 -revcomp  -maxsize 600000 -evt 0.0001
done
#

meme AC.fa -dna -oc ./meme/AC -nostatus -time 1440000 -mod zoops  -minw 6 -maxw 50 -revcomp  -maxsize 600000 -evt 0.0001
#
#
fimo --oc ./fimo/AC --verbosity 1 --thresh 1.0E-4 prodoric.meme AC.fa
#
#
#