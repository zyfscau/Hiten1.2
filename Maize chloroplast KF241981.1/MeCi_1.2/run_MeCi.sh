for sample_id in ppr329-1
do
	cutadapt -j 16 -a AGATCGGAAGAG -A AGATCGGAAGAG -m 50 --pair-filter=both -o ${sample_id}WT5-2_1.clean.fastq.gz -p ${sample_id}WT5-2_2.clean.fastq.gz  ${sample_id}WT5-2_1.gz  ${sample_id}WT5-2_2.gz
	fastqc -o ./qc_data/ -t 16 ${sample_id}_1.clean.fastq ${sample_id}_2.clean.fastq
	/mnt/e/ubuntu/MeCi-v1.0.1/bin/flash ${sample_id}_1.clean.fastq ${sample_id}_2.clean.fastq -t 28 -M 150 -o ${sample_id}
	seqkit seq -p -r ${sample_id}.notCombined_2.fastq > ${sample_id}.notCombined_2.pr.fastq
	cat ${sample_id}.extendedFrags.fastq ${sample_id}.notCombined_1.fastq ${sample_id}.notCombined_2.pr.fastq > ${sample_id}.clean.bf.fastq
	rm ${sample_id}.extendedFrags.fastq ${sample_id}.notCombined_1.fastq ${sample_id}.notCombined_2.pr.fastq ${sample_id}.notCombined_2.fastq
	cutadapt -j 16 -m 25 -o ${sample_id}.clean.fastq ${sample_id}.clean.bf.fastq
	rm ${sample_id}.clean.bf.fastq ${sample_id}_1.clean.fastq ${sample_id}_2.clean.fastq
	rm ${sample_id}.hist ${sample_id}.histogram
	perl /mnt/e/ubuntu/MeCi-v1.0.1/run_MeCi.pl --in1 ${sample_id}.clean.fastq --genome NC_007982.1.fasta --chunksize 200000 --max_process 16, --paired
#	rm ${sample_id}.clean.fastq
	rm -r ./out/01.index
	rm -r ./out/02.reads
	rm -r ./out/03.split
	rm -r ./out/04.blast
	rm -r ./out/05.parse
	mv ./out ./${sample_id}_out
done
 fastqc -f fastq -noextract -t 16 ppr329-1_1.clean.fastq.gz ppr329-1_2.clean.fastq.gz ppr329-2_1.clean.fastq.gz ppr329-2_2.clean.fastq.gz WT3-1_1.clean.fastq.gz WT3-1_2.clean.fastq.gz WT3-2_1.clean.fastq.gz WT3-2_2.clean.fastq.gz