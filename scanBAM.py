from pysam import *
import sys
from pyfaidx import Fasta

def has_mismatch_in_interval(reference, bamfile, chrom, start, end):
    """
    Return whether there is a mismatch in the interval (start, end) in any read mapping to the given chromosome.

    reference -- an IndexedFasta object
    """
    for column in bamfile.pileup(chrom, start, end):
        refbase = reference[chrom][column.pos:column.pos+1] #.decode()
        for piledup in column.pileups:
            if piledup.indel != 0:  # Insertion is positive; deletion is negative
                continue
            querybase = piledup.alignment.query_sequence[piledup.query_position]
            if refbase != querybase:
                # Mismatch
                return True
    return False

def reads_with_mismatch_in_interval(reference, bamfile, chrom, start, end):
    """
    Return reads with at least one mismatch in the interval (start, end) on the given chromosome.

    reference -- an IndexedFasta object
    """
    mreads = []
    for column in bamfile.pileup(chrom, start, end):
        refbase = reference[chrom][column.pos:column.pos+1] #.decode()
        for piledup in column.pileups:
            if piledup.indel != 0:  # Insertion is positive; deletion is negative
                continue
            querybase = piledup.alignment.query_sequence[piledup.query_position]
            if refbase != querybase:
                # Mismatch
                mreads.append(piledup.alignment.query_name)
    return mreads

def aln_has_mismatch_in_interval(reference, bamformat, aln, chrom, start, end):
    """
    Check if alignment at least one mismatch in the interval (start, end) on the given chromosome.

    reference -- an IndexedFasta object
    """
    if bamformat == 'UCSC':
        chrom = chrom[3:]
        if chrom == 'M': chrom = 'MT'

    if (aln.get_overlap(start, end)) == 0: 
        return False

    ds = start - aln.reference_start - 1
    de = aln.reference_end - end 
    qseq = aln.query_alignment_sequence[(aln.query_alignment_start+ds):(aln.query_alignment_end-de)] 
    rseq = reference[chrom][(aln.reference_start+ds):(aln.reference_end-de)] 

    if len(qseq) != len(rseq): return True
    
    if aln.qname == "SRR1576181.6093138": 
        sys.stderr.write(str(qseq)+'\n')
        sys.stderr.write(str(rseq)+'\n')

    for i in range(0, len(qseq)):
        if str(qseq[i]) != str(rseq[i]): 
            return True
    return False

# Determine if an alignment has few enough mismatches to pass.
def mismatches_ok(aln, max_mismatches=1):
    try:
        nm = aln.get_tag('nM') 
    except:
        try:
            nm = aln.get_tag('NM') 
        except:
            return(-1) # Could not find a tag for number of mismatches
    return (nm <= max_mismatches)

# Is pairing OK? Single-end passes automatically; paired-end passes if properly paired.
def pairing_ok(aln):
    if not aln.is_paired: 
        return True
    elif aln.is_proper_pair:
        return True
    else:
        return False
    
# Is the read multi-mapping? Fail if so
def multimapping_ok(aln, max_loci=1):
    try:
        if aln.get_tag('NH') > max_loci: 
            return False
        else:
            return True
    except:
        try:
            if aln.get_tag('XT') == 'R': # XT:A:U means unique, XT:A:R means repeat, apparently 
                return False 
            else: 
                return True
        except:
            return(-1) # Could not find a tag for number of loci

# Check for spliced alignments [Not used currently]
def cigar_ok(aln): 
    if 'N' in aln.cigarstring:
        return False # Discard if spliced alignment
    return True

# Get the MD:Z tag referring to mismatches [not used currently]
def get_mismatch_loc(aln):
    try: 
        mmloc = x.get_tag('MD:Z')
    except:
        return(-1)
    return mmloc
    
# Identify alignments that pass all the criteria defined above
def find_nice_alns(region, bamfile, max_mismatches=1):
    good_alns = []
    try:
        iter = bamfile.fetch(region[0],region[1],region[2])
    except:
        sys.exit("Region" + region[0] + ' ' + str(region[1]) + ' ' + str(region[2]) + '\nBAM file' + bamfile + '\nMake sure that you have an indexed BAM file!')
    for x in iter:
        if mismatches_ok(x) and pairing_ok(x) and multimapping_ok(x):        
            good_alns.append(x)
    return(good_alns)
        
# Check input arguments
if len(sys.argv) < 4:
    sys.exit("python analyse_sam.py <peptide file> <indexed FASTA file> <BAM files>")

aln_table = {} # A dictionary that will contain, for each BAM file, a dictionary of peptide alignment counts.
debug_file = open("debug.log","w")

# Parse the TSV file to keep some crucial info in memory
peptides = []
tsv_info = {'score':{}, 'psmcount':{}, 'txtype':{}, 'chrom':{},'start':{},'end':{}, 'pepcoord':{}}
with open(sys.argv[1]) as tsv:
    tsv.readline() # Skip header
    for line in tsv:
        [pepseq, pepcoord, annotation, chrom, strand, start, stop, msgf_score, psm_count] = line.rsplit()
        peptides.append(pepseq)
        tsv_info['score'][pepseq] = msgf_score
        tsv_info['psmcount'][pepseq] = psm_count
        if annotation[0:3] == 'lnc':
            tsv_info['txtype'][pepseq]='lnc'
        elif annotation[0:3] == 'PGO':
            tsv_info['txtype'][pepseq]='pg'
        else: 
            tsv_info['txtype'][pepseq]='nov'
        tsv_info['chrom'][pepseq] = chrom
        tsv_info['start'][pepseq] = start        
        tsv_info['end'][pepseq] = stop
        tsv_info['pepcoord'][pepseq] = pepcoord

# Read indexed FASTA reference file.
refFasta = Fasta(sys.argv[2])

# Go through all the BAM files. 
for bam in sys.argv[3:]:
    sys.stderr.write(bam + '\n')
    aln_count = {} # A dictionary that will collect 'good alignment' counts by peptide.
    bamfile = AlignmentFile(bam,"rb") 
    suspected_bamchromformat = 'ENS' # 'ENS' or 'UCSC' 
    for ref in bamfile.references:
        if ref.startswith('chr'): suspected_bamchromformat = 'UCSC'
    max_mismatches = 1
    for p in peptides:
        
        chrom = tsv_info['chrom'][p]
        # BAM file in UCSC format, TSV file in Ensembl
        if suspected_bamchromformat=="UCSC":
            if not chrom.startswith("chr"): chrom = 'chr' + chrom 
            if chrom == "chrMT": chrom = "chrM"
        # CASE 2: BAM file in Ensembl format, TSV file in UCSC
        else:
            if chrom.startswith("chr"): chrom = chrom[3:]
            if chrom == "M": chrom = "MT"
        
        # First check for alignments spanning the actual peptide locus. (example of coordinate format: chr6_29894499_29894540_+)
        pcoord = tsv_info['pepcoord'][p]
        (pchr, pstart, pstop, pstrand) = pcoord.split('_') # We only need the start and stop; the chromosome info has already been processed above

        # Obtain a list of reads that have a mismatch in the peptide locus
        # pep_mism_reads = reads_with_mismatch_in_interval(refFasta, bamfile, chrom, int(pstart), int(pstop))
        # print("Peptide coords: ", chrom, pstart, pstop)
        # Then find the intervals that we should check for alignments along the whole annotated locus. This could be a single interval (for a pseudogene or novel transcribed region)
        # or a set of regions (for spliced long non-coding RNAs).
        start = tsv_info['start'][p]
        stop = tsv_info['end'][p]
        regions_to_scan = []
        for i in range(0, len(start.rsplit(";"))):
            regions_to_scan.append( (chrom, int(start.rsplit(";")[i]), int(stop.rsplit(";")[i])) )        
        passed_alns = set()
        n_failed_due_to_mismatch = 0
        for r in regions_to_scan:
            good_alns = find_nice_alns(r, bamfile)
            for a in good_alns:
                if not aln_has_mismatch_in_interval(refFasta, suspected_bamchromformat, a, chrom, int(pstart), int(pstop)):
                    passed_alns.add(a)
                else:
                    debug_file.write('Peptide mismatch ' + a.qname + ' ' + pstart + ' ' + pstop + '\n')
                    n_failed_due_to_mismatch += 1
        #if not a.qname in pep_mism_reads:
                    #passed_alns.add(a.qname)
        aln_count[p] = len(passed_alns)
        debug_file.write('BAM file: ' + bam + ' Peptide locus ' + chrom + ':' + pstart + '-' + pstop + '\n')
        debug_file.write('Passed through ' + str(len(passed_alns)) + ' alignments. Found ' +  str(n_failed_due_to_mismatch) + ' alignments with mismatches in peptide locus' + '\n#####################\n')
    aln_table[bam]=aln_count

bam = sorted(list(aln_table.keys()))

# Write output file header.
sys.stdout.write('sequence\ttxtype\tscore\tpsmcount\t')
for i in range(0,len(bam)):
        sys.stdout.write(bam[i].split('.')[0].split('/')[-1])
        if (i < (len(bam)-1)): sys.stdout.write('\t')
        else: sys.stdout.write('\n')

# And the counts.
for pep in sorted(aln_table[bam[0]].keys()):
    sys.stdout.write(pep + '\t' + tsv_info['txtype'][pep] + '\t' + tsv_info['score'][pep] + '\t' + tsv_info['psmcount'][pep] + '\t')
    for i in range(0,len(bam)): 
        sys.stdout.write(str(aln_table[bam[i]][pep]))
        if (i < (len(bam)-1)): sys.stdout.write('\t')
    sys.stdout.write('\n')


