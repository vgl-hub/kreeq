validate -f testFiles/to_correct.fasta -r testFiles/to_correct.fastq -o vcf
embedded
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
sequence2	25	.	a	T	42	PASS	.	GT	1/1
sequence2	65	.	t	C	42	PASS	.	GT	1/1
sequence3	33	.	Ga	G,GCa	22	PASS	.	GT	1/1
sequence4	33	.	GA	GCA	40	PASS	.	GT	1/1
sequence5	69	.	Aa	A	36	PASS	.	GT	1/1
sequence6	68	.	AT	AAT	36	PASS	.	GT	1/1
sequence7	25	.	a	T	16	PASS	.	GT	1/1
sequence7	35	.	t	A	42	PASS	.	GT	1/1
sequence8	72	.	c	A	16	PASS	.	GT	1/1
sequence9	50	.	t	C,T	2	PASS	.	GT	1/1
sequence10	25	.	a	T	0	PASS	.	GT	1/1
sequence10	27	.	c	G	42	PASS	.	GT	1/1
sequence11	25	.	a	T	0	PASS	.	GT	1/1
sequence11	27	.	c	G	42	PASS	.	GT	1/1
sequence11	65	.	t	C	42	PASS	.	GT	1/1
sequence12	25	.	a	T	0	PASS	.	GT	1/1
sequence12	27	.	c	G	12	PASS	.	GT	1/1
sequence12	35	.	t	A	42	PASS	.	GT	1/1
sequence13	25	.	a	T	42	PASS	.	GT	1/1
sequence13	68	.	AT	AAT	36	PASS	.	GT	1/1
sequence14	33	.	GA	GCA	40	PASS	.	GT	1/1
sequence14	67	.	AT	AAT	36	PASS	.	GT	1/1
sequence15	46	.	AT	AAT	36	PASS	.	GT	1/1
sequence16	66	.	GA	GAA	-15	PASS	.	GT	1/1
sequence17	25	.	a	T	42	PASS	.	GT	1/1
sequence17	69	.	Aa	A	36	PASS	.	GT	1/1
sequence18	33	.	Ga	G,GCa	22	PASS	.	GT	1/1
sequence18	70	.	Aa	A	36	PASS	.	GT	1/1
sequence19	33	.	Ga	G,GCa	10	PASS	.	GT	1/1
sequence19	43	.	Ta	T	40	PASS	.	GT	1/1
sequence20	68	.	Aa	A	-15	PASS	.	GT	1/1
sequence21	24	.	Ca	C	2	PASS	.	GT	1/1
