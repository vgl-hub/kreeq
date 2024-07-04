kreeq validate -f testFiles/to_correct.fasta -r testFiles/to_correct.fastq -o vcf --search-depth 50 --max-span 32
embedded
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
sequence2	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence2	65	.	t	C	0	PASS	.	GT:GQ	1/1:0
sequence3	33	.	Ga	G	0	PASS	.	GT:GQ	1/1:0
sequence4	33	.	GA	GCA	0	PASS	.	GT:GQ	1/1:0
sequence5	69	.	Aa	A	0	PASS	.	GT:GQ	1/1:0
sequence6	68	.	AT	AAT	0	PASS	.	GT:GQ	1/1:0
sequence7	25	.	aCGTACATGCt	TCGTACATGCA	0	PASS	.	GT:GQ	1/1:0
sequence8	72	.	c	A	0	PASS	.	GT:GQ	1/1:0
sequence9	50	.	t	C	0	PASS	.	GT:GQ	1/1:0
sequence10	25	.	aCc	TCG	0	PASS	.	GT:GQ	1/1:0
sequence11	25	.	aCc	TCG	0	PASS	.	GT:GQ	1/1:0
sequence11	65	.	t	C	0	PASS	.	GT:GQ	1/1:0
sequence12	25	.	aCcTACATGCt	TCGTACATGCA	0	PASS	.	GT:GQ	1/1:0
sequence13	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence13	68	.	AT	AAT	0	PASS	.	GT:GQ	1/1:0
sequence14	33	.	GA	GCA	0	PASS	.	GT:GQ	1/1:0
sequence14	67	.	AT	AAT	0	PASS	.	GT:GQ	1/1:0
sequence15	46	.	AT	AAT	0	PASS	.	GT:GQ	1/1:0
sequence15	67	.	AT	AAT	0	PASS	.	GT:GQ	1/1:0
sequence16	67	.	AT	AAAT	0	PASS	.	GT:GQ	1/1:0
sequence17	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence17	69	.	Aa	A	0	PASS	.	GT:GQ	1/1:0
sequence18	33	.	Ga	G	0	PASS	.	GT:GQ	1/1:0
sequence18	70	.	Aa	A	0	PASS	.	GT:GQ	1/1:0
sequence19	34	.	aCAGTGATGTa	TGCAGTGATGT	0	PASS	.	GT:GQ	1/1:0
sequence20	69	.	Aa	A	0	PASS	.	GT:GQ	1/1:0
sequence21	25	.	at	TC	0	PASS	.	GT:GQ	1/1:0
sequence22	30	.	g	C	0	PASS	.	GT:GQ	1/1:0
sequence22	70	.	g	T	0	PASS	.	GT:GQ	1/1:0
sequence23	40	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence23	75	.	c	G	0	PASS	.	GT:GQ	1/1:0
