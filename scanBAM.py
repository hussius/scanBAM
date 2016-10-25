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

def parse_coords(cstr):
    return(cstr.split('_')[:-1])

def is_in_interval(coord, interval_start, interval_end):
    return ( ( int(coord) >= int(interval_start)) and (int(coord) <= int(interval_end) ) )

def get_overlap(s1, e1, s2, e2):
    """ 
    Get the coordinates of the overlap between two intervals
    """
    if s1 > e2 or e1 < s2: return None
    if s1 <= s2 and e1 <= e2: return (s2, e1) # Will also work for s1 == s2 and e1 == e2 
    if s1 <= s2 and e1 >= e2: return (s2, e2) # Alignment enclosed in peptide
    if s1 >= s2 and e1 <= e2: return (s1, e1) # Peptide enclosed in alignment
    if s1 >= s2 and e1 >= e2: return (s1, e2)
    sys.exit('Check your numbers')

def get_peptide_segments(pep_info, p):
    """
    In some cases, peptides come from spliced segments. In these cases we need to 
    infer the peptide coding regions from the tsv file.
    """
    start = pep_info['start'][p]
    end = pep_info['end'][p]
    chr = pep_info['chrom'][p]
    pepcoor = pep_info['pepcoord'][p]
    pep_coord = parse_coords(pepcoor)    
    if len(start.split(';')) == 1: # Not spliced
        return([pep_coord[1], pep_coord[2]])
    # If spliced
    estarts = start.split(';')
    eends = end.split(';')
    #print(chr)
    #print(pep_coord)
    assert(pep_coord[0]==chr)
    pep_start_coord_exon = -1
    pep_end_coord_exon = -1
    for i in range(0, len(estarts)):
        assert(int(estarts[i]) < int(eends[i]))
        if is_in_interval(pep_coord[1], estarts[i], eends[i]): pep_start_coord_exon = i
        if is_in_interval(pep_coord[2], estarts[i], eends[i]): pep_end_coord_exon = i
        # For peptides where start or end cannot be assigned to an exon, assume it belongs to the same one as the end that can be assigned.
        # This is a rare case
    #print('Peptide coordinates:', pep_coord)
    #print('Inferred start exon: ' + str(pep_start_coord_exon) + ' (' + str(estarts[pep_start_coord_exon]) + '-' + str(eends[pep_start_coord_exon]))
    #print('Inferred end exon: ' + str(pep_end_coord_exon) + ' (' + str(estarts[pep_end_coord_exon]) + '-' + str(eends[pep_end_coord_exon]))
    #input()
    if pep_start_coord_exon == -1 or pep_end_coord_exon == -1:
        if pep_start_coord_exon == -1: pep_start_coord_exon = pep_end_coord_exon
        if pep_end_coord_exon == -1: pep_end_coord_exon = pep_start_coord_exon
    if pep_start_coord_exon == pep_end_coord_exon:
        # Peptide contained within one exon. Can print out the original coordinates.
        return([pep_coord[1], pep_coord[2]])
    else:
        # Need to stitch together the regions from annotation.
        # First entry, for the starting exon
        # Chromosome, peptide start coord, end of the starting exon, etc
        entries = []
        entries.append( [pep_coord[1], eends[pep_start_coord_exon]] )
        # Second entry, for the ending exon
        entries.append( [estarts[pep_end_coord_exon], pep_coord[2]])
        return(entries)

def aln_has_mismatch_in_interval(reference, bamformat, aln, chrom, start, end):
    """
    Check if (unspliced) alignment has at least one mismatch in the interval (start, end) on the given chromosome.
    
    reference -- an IndexedFasta object
    """

    qseq = ''
    rseq = ''

    dbg = False
    #if 'SRR1027172.10275983' in aln.query_name: dbg = True
    #if dbg: print('### DEBUG ###\t' + str(start) + '\t' + str(end) + '\t' + aln.query_name)

    if bamformat == 'UCSC':
        chrom = chrom[3:]
        if chrom == 'M': chrom = 'MT'

    if (aln.get_overlap(start, end)) == 0: 
        return False

    #if dbg:
    #    print('Alignment reference start:\t' + str(aln.reference_start))
    #    print('Alignment reference end:\t' + str(aln.reference_end))

    if start >= aln.reference_start and end <= aln.reference_end:
        # CASE 1: Peptide contained within aligned segment
        # Peptide                            start ------ end
        # Alignment    aln.reference_start ----------------------- aln.reference_end
        ds = start - aln.reference_start # Offset on aligned sequence, from its start
        de = aln.reference_end - end     # Offest on aligned sequence, from its end
        qseq = aln.query_alignment_sequence[ds:(aln.query_alignment_end-aln.query_alignment_start-de)] 
        rseq = reference[chrom][(aln.reference_start+ds):(aln.reference_end-de)] 
    
    elif start < aln.reference_start:
        # Peptide starts before aligned segment
        if end <= aln.reference_end:
            # CASE 2: Peptide left-overlapping aligned segment
            # Peptide              start ---------- end
            # Alignment  aln.reference_start ---------------- aln.reference_end
            ds = end - aln.reference_start 
            # Offset on aligned sequence, from the start
            qseq = aln.query_alignment_sequence[:ds] 
            rseq = reference[chrom][aln.reference_start:(aln.reference_start+ds)]
        else: 
            # CASE 3: (rare) alignment contained within peptide. start < ref start and end > ref end
            # Peptide            start --------------------------------- end
            # Alignment aln.reference_start --------------------- aln.reference_end
            qseq = aln.query_alignment_sequence
            rseq = reference[chrom][aln.reference_start:aln.reference_end] 
    else: # start >= ref start and end >= ref end
        # CASE 4: Peptide right-overlapping aligned segment
        # Peptide                                    start ---------------- end
        # Alignment     aln.reference_start --------------------- aln.reference_end
        assert start >= aln.reference_start and end > aln.reference_end
        de = aln.reference_end - start 
        qseq = aln.query_alignment_sequence[(aln.query_alignment_end-aln.query_alignment_start-de):(aln.query_alignment_end-aln.query_alignment_start)] 
        rseq = reference[chrom][(aln.reference_end-de):(aln.reference_end)] 
    
    if qseq == '' or rseq == '':
        #print(aln)
        sys.exit('Could not extract sequence.')

    if len(qseq) != len(rseq): 
        # Should happen for insertions and deletions only
        debug_file.write('Indel: ' + aln.cigarstring + '\n')
        return True
    
    for i in range(0, len(qseq)):
        if str(qseq[i]) != str(rseq[i]):
            debug_file.write('Mismatch: ' + '\n' + str(qseq) + '\n' + str(rseq) + '\n')
            # More details
            debug_file.write("Start of peptide: " + '\t' + str(start) + '\n')
            debug_file.write("End of peptide: " + '\t' + str(end) + '\n')
            debug_file.write("Start of alignment: " + '\t' + str(aln.reference_start) + '\n')
            debug_file.write("End of alignment: " + '\t' + str(aln.reference_end) + '\n')
            debug_file.write(str(aln))

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

# Check for spliced alignments
# UNDER DEVELOPMENT :)
def cigar_ok(aln): 
    if 'N' in aln.cigarstring:
        # Find out where the read hits the reference
        ct = aln.cigartuples
        curr_loc = aln.reference_start
        aln_starts = []
        aln_ends = []
        for tup in ct:
            if tup[0] == 0:
                aln_starts.append(curr_loc)
                aln_ends.append(curr_loc + tup[1])
            # This needs to be updated in any case    
            curr_loc += tup[1]
        for e in zip(aln_starts, aln_ends):
            mm = aln_has_mismatch_in_interval(reference, bamformat, aln, chrom, e[0], e[1])
        input()
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
        
def find_fusion_alns(regions, bamfile):
    """
    Try to find alignments corresponding to a fusion event.
    Right now just deals with 2 regions.
    """ 
    fus_alns = []
    chrom1, start1, end1 = regions[0]
    try:
        chrom2, start2, end2 = regions[1]
    except:
        print("Abnormal exit")
        print(regions)
        sys.exit(0)
    #print('Partner should be at:' + str(chrom2) + ':' + str(start2) + '-' + str(end2))
    iter = bamfile.fetch(chrom1, start1, end1)
    for x in iter:
        if bamfile.getrname(x.reference_id) == chrom1 and bamfile.getrname(x.next_reference_id) == chrom2:
            if (x.next_reference_start >= start2 and x.next_reference_start <= end2):
                if mismatches_ok(x) and pairing_ok(x) and multimapping_ok(x):
                    fus_alns.append(x)
    return(fus_alns)

def compare_seqs(aln, pcoords, reference, chrformat, chrom):

    if chrformat == 'UCSC':
        chrom = chrom[3:]
        if chrom == 'M': chrom = 'MT'    
    mm = 0
    # Step through alignment entries
    ct = aln.cigartuples
    curr_loc = aln.reference_start
    offset_in_read = 0
    for tup in ct:
        if tup[0] == 0: # Match
            aln_seg_start = curr_loc
            aln_seg_end = curr_loc + tup[1]
            for seg in pcoords:
                ol = get_overlap(int(seg[0]),int(seg[1]),aln_seg_start,aln_seg_end)
                if ol:
                    overlap_length = ol[1] - ol[0] + 1
                    overlap_offset = ol[0] - curr_loc
                    qseq = aln.query_sequence[(aln.query_alignment_start+offset_in_read+overlap_offset):(aln.query_alignment_start+offset_in_read+overlap_offset+overlap_length)]
                    rseq = reference[chrom][ol[0]:(ol[1]+1)]
                    assert( overlap_offset >= 0)
                    for i in range(0, len(qseq)): # It can happen that len(rseq) > len(qseq) if we are at the end of the read, but that's ok! We are only interested in mismatches in qseq w r t rseq
                        if str(qseq[i]) != str(rseq[i]):
                            mm += 1
                            debug_file.write('Mismatch in spliced segment: ' + '\n' + str(qseq) + '\n' + str(rseq) + '\n')
            offset_in_read += tup[1] # Keep track of location in read
        curr_loc += tup[1] # Keep track of location on reference
    return(mm)

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
        elif annotation[0:3] == "Fus":
            tsv_info['txtype'][pepseq]='fus'
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

        dbg = False
        if p == "AAAEIDEEPVSK": dbg = True
        chroms = tsv_info['chrom'][p].split(';')
        # Coordinates for the actual peptide locus or loci (in the case of spliced peptides). (example of coordinate format: chr6_29894499_29894540_+)
        pcoords = get_peptide_segments(tsv_info, p)
        # Below commented out when implementing spliced peptide handling
        #pcoord = tsv_info['pepcoord'][p]
        #(pchr, pstart, pstop, pstrand) = pcoord.split('_') # We only need the start and stop; the chromosome info has already been processed above
        # Find the intervals that we should check for alignments along the whole annotated locus. This could be a single interval (for a pseudogene or novel transcribed region)
        # or a set of regions (for spliced long non-coding RNAs).
        start = tsv_info['start'][p]
        stop = tsv_info['end'][p]
        regions_to_scan = []
        for i in range(0, len(start.rsplit(";"))):
            if len(chroms) > 1:
                chrom = chroms[i]
            else:
                chrom = chroms[0]
            # CASE 1: BAM file in UCSC format, TSV file in Ensembl
            if suspected_bamchromformat=="UCSC":
                if not chrom.startswith("chr"): chrom = 'chr' + chrom 
                if chrom == "chrMT": chrom = "chrM"
            # CASE 2: BAM file in Ensembl format, TSV file in UCSC
            else:
                if chrom.startswith("chr"): chrom = chrom[3:]
                if chrom == "M": chrom = "MT"
            regions_to_scan.append( (chrom, int(start.rsplit(";")[i]), int(stop.rsplit(";")[i])) )        
        passed_alns = set()
        n_failed_due_to_mismatch = 0

        ###
        # Right now we just search for evidence separately in the putative fusion regions
        # I e not requiring a chimeric alignment
        # So the code below is commented out
        ###

        #if tsv_info['txtype'][p] == 'fus' and len(regions_to_scan) > 1:
        #    fus_alns = []
        #    assert(len(regions_to_scan)%2==0)
        #    for i in range(int(len(regions_to_scan)/2)):
        #        fus_alns.append( find_fusion_alns(regions_to_scan[2*i:(2*i+2)], bamfile))
        for r in regions_to_scan:
            print(r)
            input()
            good_alns = find_nice_alns(r, bamfile)
            for a in good_alns:
                # Check if spliced alignment
                if 'N' in a.cigarstring:
                    # If spliced, find out the aligned bits
                    #print("Spliced alignment")
                    #print(str(a.reference_start) + '\t' + a.cigarstring)
                    ct = a.cigartuples
                    curr_loc = a.reference_start
                    aln_starts = []
                    aln_ends = []
                    for tup in ct:
                        if tup[0] == 0:
                            aln_starts.append(curr_loc)
                            aln_ends.append(curr_loc + tup[1]-1)
                        # This needs to be updated in any case    
                        curr_loc += tup[1]-1
                    overlap = False
                    #if dbg: print('Spliced')
                    for e in zip(aln_starts, aln_ends):
                        for seg in pcoords:
                            ol = get_overlap(int(seg[0]), int(seg[1]), e[0], e[1])
                            if ol:
                                overlap = True
                                #print('Overlap between: ' + str(seg) + ' and ' + str(e) + ': ' + str(ol))
                    if overlap: 
                        # Check if they have mismatches
                        mm = compare_seqs(a, pcoords, refFasta, suspected_bamchromformat, chrom) 
                        if mm == 0:
                            passed_alns.add(a)
                # Not spliced alignment - just make sure there is no mismatch in the peptide region
                else:
                    mm = False 
                    for seg in pcoords: 
                        if aln_has_mismatch_in_interval(refFasta, suspected_bamchromformat, a, chrom, int(seg[0]), int(seg[1])):
                            mm = True
                    if not mm:
                        passed_alns.add(a)
                    else:
                        debug_file.write('Peptide mismatch ' + a.qname + ' ' + seg[0] + ' ' + seg[1] + '\n')
                        n_failed_due_to_mismatch += 1
                    
        #if not a.qname in pep_mism_reads:
                    #passed_alns.add(a.qname)
        aln_count[p] = len(passed_alns)
        
        debug_file.write('BAM file: ' + bam + ' Peptide locus ' + chrom + ':' + pcoords[0][0] + '-' + pcoords[0][1] + '\n')
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



