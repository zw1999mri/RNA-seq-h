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


ls 3-bam/*.bam |while read id;do (featureCounts -T 30 -t gene -g gene_id -p --countReadPairs -a *.gtf -o 4-count/$(basename ${id} ".bam").txt ${id});done








