# Troubleshooting

## TrimGalore Error
Check the log file. It's most likely because that fastq has no reads. <br/>
If that is the case, one quick fix is to remove empty fastqs from the SAMPLES list defined in the Snakefile. <br/>
One way to find the fastqs with no reads is to run:
```
find . -type f -name "*.gz" -size 20c
```
in the directory containing your demultiplexed fastq.gz files.

Once you have a newline delimited list, you can add this block to the Snakefile:
```
empty_fastqs = open('empty_fastq_list', 'r')
for fastq in empty_fastqs.readlines():
	SAMPLES.remove(fastq.strip())
empty_fastqs.close()
```
It is important that you add this just after the block where SAMPLES is created and filled with sample names from the sample sheet. <br/>
The pipeline should then only process the fastqs which aren't empty.


**Note that this is __NOT__ an ideal solution and the empty fastqs will not be represented in any of the summary stat files.**
