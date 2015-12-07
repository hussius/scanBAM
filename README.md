# scanBAM
This is a script that searches a set of BAM files for RNA-seq alignments that match a provided set of peptide coordinates.

The intended use case is when a set of "unusual" peptides (i.e., peptides that do not correspond to annotated coding regions in the genome) has been detected by high-throughput mass spectrometry and we want to check public RNA-seq data sets for evidence of the corresponding genomic regions being transcribed. 

<h3>Dependencies</h3>
- pysam
- pyfaidx 

<h3>ToDo</h3>
Add support for suspected fusion peptides
