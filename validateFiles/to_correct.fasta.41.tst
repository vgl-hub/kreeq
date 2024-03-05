validate -f testFiles/to_correct.fasta -r testFiles/to_correct.fastq -o vcf
embedded
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
sequence2	25	.	a	T	42	PASS	.	GT:GQ	1/1:42
sequence2	65	.	t	C	42	PASS	.	GT:GQ	1/1:42
sequence3	33	.	Ga	G,GCa	22	PASS	.	GT:GQ	1/1:42
sequence4	33	.	GA	GCA	40	PASS	.	GT:GQ	1/1:40
sequence5	69	.	Aa	A	36	PASS	.	GT:GQ	1/1:36
sequence6	68	.	AT	AAT	36	PASS	.	GT:GQ	1/1:36
sequence7	25	.	a	T	16	PASS	.	GT:GQ	1/1:16
sequence7	35	.	t	A	42	PASS	.	GT:GQ	1/1:42
sequence8	72	.	c	A	16	PASS	.	GT:GQ	1/1:16
sequence9	50	.	t	C,T	19	PASS	.	GT:GQ	1/1:38
sequence10	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence10	27	.	c	G	42	PASS	.	GT:GQ	1/1:42
sequence11	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence11	27	.	c	G	42	PASS	.	GT:GQ	1/1:42
sequence11	65	.	t	C	42	PASS	.	GT:GQ	1/1:42
sequence12	25	.	a	T	0	PASS	.	GT:GQ	1/1:0
sequence12	27	.	c	G	12	PASS	.	GT:GQ	1/1:12
sequence12	35	.	t	A	42	PASS	.	GT:GQ	1/1:42
sequence13	25	.	a	T	42	PASS	.	GT:GQ	1/1:42
sequence13	68	.	AT	AAT	36	PASS	.	GT:GQ	1/1:36
sequence14	33	.	GA	GCA	40	PASS	.	GT:GQ	1/1:40
sequence14	67	.	AT	AAT	36	PASS	.	GT:GQ	1/1:36
sequence15	46	.	AT	AAT	36	PASS	.	GT:GQ	1/1:36
sequence16	66	.	GA	GAA	2	PASS	.	GT:GQ	1/1:2
sequence17	25	.	a	T	42	PASS	.	GT:GQ	1/1:42
sequence17	69	.	Aa	A	36	PASS	.	GT:GQ	1/1:36
sequence18	33	.	Ga	G,GCa	22	PASS	.	GT:GQ	1/1:42
sequence18	70	.	Aa	A	36	PASS	.	GT:GQ	1/1:36
sequence19	33	.	Ga	G,GCa	10	PASS	.	GT:GQ	1/1:20
sequence19	43	.	Ta	T	40	PASS	.	GT:GQ	1/1:40
sequence20	68	.	Aa	A	2	PASS	.	GT:GQ	1/1:2
sequence21	24	.	Ca	C	2	PASS	.	GT:GQ	1/1:2
