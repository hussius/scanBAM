# scanBAM
This is a script that searches a set of BAM files for RNA-seq alignments that match a provided set of peptide coordinates. "Match" here means that an alignment is spanning the peptide region, that it does not have any mismatches in the peptide region itself (although it can have mismatches elsewhere), that it is not multi-mapping, and that it is properly paired (in the paired-end case).

The intended use case is when a set of "unusual" peptides (i.e., peptides that do not correspond to annotated coding regions in the genome) has been detected by high-throughput mass spectrometry and we want to check public RNA-seq data sets for evidence of the corresponding genomic regions being transcribed. 

I might add some example input files here if there is any interest.

<h3>Dependencies</h3>
- pysam
- pyfaidx 

<h3>ToDo</h3>
Add support for suspected fusion peptides

Add example input files
