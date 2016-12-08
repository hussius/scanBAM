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

def get_peptide_segments(pep_info, p, suspected_format):
    """
    In some cases, peptides come from spliced segments. In these cases we need to 
    infer the peptide coding regions from the tsv file.
    """
    start = pep_info['start'][p]
    end = pep_info['end'][p]
    chr = pep_info['chrom'][p]
    pepcoor = pep_info['pepcoord'][p]
    pep_coord = parse_coords(pepcoor)    
    #print(pep_info['chrom'][p])
    #print(pepcoor)
    assert(pep_coord[0]==chr)
    pep_chr = pep_coord[0]
    # CASE 1: BAM file in UCSC, TSV in ENSEMBL
    if suspected_format=="UCSC":
        if not pep_chr.startswith("chr"): pep_chr = 'chr' + pep_chr
        if pep_chr == "chrMT": pep_chr = "chrM"
    # CASE 2: BAM file in Ensembl format, TSV file in UCSC
    else:
        if pep_chr.startswith("chr"): pep_chr = pep_chr[3:]
        if pep_chr == "M": pep_chr = "MT"

    if len(start.split(';')) == 1: # Not spliced
        return([(pep_chr, int(pep_coord[1]), int(pep_coord[2]) )])
    # If spliced
    estarts = start.split(';')
    eends = end.split(';')
    #print(chr)
    #print(pep_coord)
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
        return([(pep_chr, int(pep_coord[1]), int(pep_coord[2]) )])
    else:
        # Need to stitch together the regions from annotation.
        # First entry, for the starting exon
        # Chromosome, peptide start coord, end of the starting exon, etc
        entries = []
        entries.append( (pep_chr, int(pep_coord[1]), int(eends[pep_start_coord_exon])) )
        # Second entry, for the ending exon
        entries.append( (pep_chr, int(estarts[pep_end_coord_exon]), int(pep_coord[2])) )
        return(entries)

def aln_has_mismatch_in_interval(reference, bamformat, aln, chrom, start, end):
    """
    Check if an (unspliced) alignment has at least one mismatch in the interval (start, end) on the given chromosome.
    
    reference -- an IndexedFasta object
    """

    qseq = ''
    rseq = ''

    dbg = False

    if bamformat=="UCSC": # Need to give chromosome name that fits with ENSEMBL for now
        chrom = chrom[3:]
        if chrom == "M": chrom = "MT"

    if (end-start) > 1000:
        print('Long peptide')
        print(str(start) + ' ' + str(end) + ' ' + str(end-start))
        print('Reference sequence of peptide:\t' + str(reference[chrom][start:end]))
        input()

    if (aln.get_overlap(start, end)) == 0: 
        return False

    if dbg:
        pass
        #print('Alignment reference start:\t' + str(aln.reference_start))
        #print('Alignment reference end:\t' + str(aln.reference_end))
        #print('Peptide chromosome:\t' + chrom)
        #print('Peptide start:\t' + str(start))
        #print('Peptide end:\t' + str(end))
        #print(aln) 
        #print('Sequence of whole alignment    :\t' + aln.query_alignment_sequence)
        #print('Reference sequence of alignment:\t' + str(reference[chrom][aln.reference_start:aln.reference_end]))
        #print('Reference sequence of peptide:\t' + str(reference[chrom][start:end]))

    if start >= aln.reference_start and end <= aln.reference_end:
        if dbg: sys.stderr.write('Case 1: peptide contained within aligned segment\n')
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
            if dbg: sys.stderr.write('Case 2: peptide left-overlapping aligned segment\n')
            # CASE 2: Peptide left-overlapping aligned segment
            # Peptide              start ---------- end
            # Alignment  aln.reference_start ---------------- aln.reference_end
            ds = end - aln.reference_start 
            # Offset on aligned sequence, from the start
            qseq = aln.query_alignment_sequence[:ds] 
            rseq = reference[chrom][aln.reference_start:(aln.reference_start+ds)]
        else: 
            if dbg: sys.stderr.write('Case 3: alignment contained within peptide\n')
            # CASE 3: (rare) alignment contained within peptide. start < ref start and end > ref end
            # Peptide            start --------------------------------- end
            # Alignment aln.reference_start --------------------- aln.reference_end
            qseq = aln.query_alignment_sequence
            rseq = reference[chrom][aln.reference_start:aln.reference_end] 
    else: # start >= ref start and end >= ref end
        # CASE 4: Peptide right-overlapping aligned segment
        # Peptide                                    start ---------------- end
        # Alignment     aln.reference_start --------------------- aln.reference_end
        if dbg: sys.stderr.write('Case 4: peptide right-overlapping aligned segment\n')
        assert start >= aln.reference_start and end > aln.reference_end
        de = aln.reference_end - start 
        qseq = aln.query_alignment_sequence[(aln.query_alignment_end-aln.query_alignment_start-de):(aln.query_alignment_end-aln.query_alignment_start)] 
        rseq = reference[chrom][(aln.reference_end-de):(aln.reference_end)] 
    if dbg: sys.stderr.write(str(qseq) + '\n')
    if dbg: sys.stderr.write(str(rseq) + '\n')

    if len(qseq) != len(rseq): 
        # Should happen for insertions and deletions only
        debug_file.write('Indel: ' + aln.cigarstring + '\n')
        return True

    if qseq == '' or rseq == '':
        sys.exit('Could not extract sequence.')
    
    for i in range(0, len(qseq)):
        if str(qseq[i]).upper() != str(rseq[i]).upper():
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

def compare_seqs(aln, seg, reference, chrformat, chrom):
    return False
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
            ol = get_overlap(int(seg[1]),int(seg[2]),aln_seg_start,aln_seg_end)
            if ol:
                print('Peptide: ' + str(seg[1]) + '-' + str(seg[2]))
                print('Aln: ' + str(aln_seg_start) + '-' + str(aln_seg_end))

                overlap_length = ol[1] - ol[0] + 1
                overlap_offset = ol[0] - curr_loc

                print('Overlap length: ' + str(overlap_length))
                print('Overlap: ')
                print(ol)
                qseq = aln.query_sequence[(aln.query_alignment_start+offset_in_read+overlap_offset):(aln.query_alignment_start+offset_in_read+overlap_offset+overlap_length)]
                rseq = reference[chrom][ol[0]:(ol[1]+1)]
                print(rseq)
                print(aln.query_sequence)
                assert( overlap_offset >= 0)
                for i in range(0, len(qseq)): # It can happen that len(rseq) > len(qseq) if we are at the end of the read, but that's ok! We are only interested in mismatches in qseq w r t rseq
                    if str(qseq[i]) != str(rseq[i]):
                        mm += 1
                        print('Mismatch in spliced segment: ' + '\n' + str(qseq) + '\n' + str(rseq) + '\n')
                        #input()
                        debug_file.write('Mismatch in spliced segment: ' + '\n' + str(qseq) + '\n' + str(rseq) + '\n')
            offset_in_read += tup[1] # Keep track of location in read
        curr_loc += tup[1] # Keep track of location on reference
    return(mm)

def mismatch_in_spliced(aln, pcoords, reference, suspected_bamchromformat):
    #print(aln)
    #print(pcoords)
    mm = 0
    chrom, pstart, pend = pcoords
    if suspected_bamchromformat == 'UCSC':
        chrom = chrom[3:]
        if chrom == 'M': chrom = 'MT'
    whole_qseq = aln.query_sequence
    #print(pcoords)
    rseq = reference[chrom][pstart:pend] 
    #print('============')
    #print('Whole query sequence: ' + str(aln.query_sequence))
    #print('Peptide-corresponding reference sequence: ' + str(rseq))
    curr_loc = aln.reference_start
    offset_in_read = 0
    ct = aln.cigartuples
    for tup in ct:
        #print('CIGAR string:', aln.cigarstring, str(tup))
        if tup[0]==0: # match
            aln_seg_start = curr_loc
            aln_seg_end = curr_loc + tup[1]
            sequence_of_seg = aln.query_sequence[offset_in_read:offset_in_read+tup[1]]
            #print('Sequence of this segment: ' + ' '*offset_in_read + sequence_of_seg )
            ol = get_overlap(pstart, pend, aln_seg_start,aln_seg_end)
            #print('Peptide coordinates: ' + str(pstart) + '-' + str(pend))
            #print('Aligned coordinates: ' + str(aln_seg_start) + '-' + str(aln_seg_end))
            if ol: 
                overlap_length = ol[1] - ol[0]
                overlap_offset_aln = ol[0] - curr_loc
                #print('Overlap: ')
                #print(ol)
                #print('Overlap offset for alignment: ' + str(overlap_offset_aln))
                #print('Overlap offset for peptide: ' + str(overlap_offset_pep))
                #print('Overlap length: ' + str(overlap_length))
                aln_bit = aln.query_sequence[offset_in_read+overlap_offset_aln:(offset_in_read+overlap_offset_aln+overlap_length)]
                try:
                    pep_rseq = reference[chrom][ol[0]:ol[1]]
                    #sys.stderr.write('aligned bit:          ' + aln_bit + '\n')
                    #sys.stderr.write('matching peptide seq: ' + str(pep_rseq) + '\n')
                    for i in range(0, len(aln_bit)):
                     if aln_bit[i].lower() != str(pep_rseq[i]).lower():
                         mm += 1
                        #print('Mismatch in spliced segment')
#                        print('aligned bit:          ' + aln_bit)
#                        print('matching peptide seq: ' + str(pep_rseq))
                        #print(ol)
                        #input()
                except:
                    pass # No overlap
            else:
                pass
                #print('No overlap')
            #input()
            offset_in_read += tup[1]
        elif tup[0]==4:
            offset_in_read += tup[1]
        if tup[0] != 4: curr_loc += tup[1]
    return mm

####
#
# Check input arguments
if len(sys.argv) < 4:
    sys.exit("python analyse_sam.py <peptide file> <indexed FASTA file> <BAM files>")

aln_table = {} # A dictionary that will contain, for each BAM file, a dictionary of peptide alignment counts.
debug_file = open("debug.log","w")

# Start by parsing the TSV file 
peptides = []
tsv_info = {'score':{}, 'psmcount':{}, 'txtype':{}, 'chrom':{},'start':{},'end':{}, 'pepcoord':{}}
with open(sys.argv[1]) as tsv:
    tsv.readline() # Skip header
    for line in tsv:
        try:
            [pepseq, pepcoord, annotation, chrom, strand, start, stop, msgf_score, psm_count] = line.rsplit('\t')
        except:
            print(line.strip())
        peptides.append(pepseq) # Peptide sequence
        tsv_info['score'][pepseq] = msgf_score # MSGF score
        tsv_info['psmcount'][pepseq] = psm_count # PSM count
        # Type of transcript annotated for peptide regions
        if annotation[0:3] == 'lnc':
            tsv_info['txtype'][pepseq]='lnc'
        elif annotation[0:3] == 'PGO':
            tsv_info['txtype'][pepseq]='pg'
        elif annotation[0:3] == "Fus":
            tsv_info['txtype'][pepseq]='fus'
        else: 
            tsv_info['txtype'][pepseq]='nov'

        tsv_info['chrom'][pepseq] = chrom # The chromosome for the annotated transcript
        tsv_info['start'][pepseq] = start # Start coordinate of annotated transcript        
        tsv_info['end'][pepseq] = stop # End coordinate of annotated transcript
        tsv_info['pepcoord'][pepseq] = pepcoord # Coordinates for the peptide, e g chr6_31620200_31620244_-

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
    max_mismatches = 1 # Refers to the maximum allowed number of mismatched *for the whole RNA alignment* (no mismatched are allowed in the peptide-overlapping region)
    for p in peptides:
        # Coordinates for the actual peptide locus or loci (in the case of spliced peptides). (example of coordinate format: chr6_29894499_29894540_+)
        # If the peptides are spliced the coordinates of the segments must be inferred from the transcript annotation
        pcoords = get_peptide_segments(tsv_info, p, suspected_bamchromformat)
        # It is for these regions that we want to find overlapping alignments.
        passed_alns = set()
        n_failed_due_to_mismatch = 0

        for r in pcoords: # For each peptide segment (usually 1)
            # Find alignments overlapping this segment that have a maximal number of mismatches, are not multimapping and paired if paired-end seq is used
            good_alns = find_nice_alns(r, bamfile)
            #sys.stderr.write('Found ' + str(len(good_alns)) + ' good alignments' + '\n')
            spliced = 0
            unspliced = 0
            not_overlapping_peptide = 0
            for a in good_alns:
                # Check if spliced alignment
                if 'N' in a.cigarstring:
                    spliced += 1
                    # If spliced, find out the aligned bits
                    ct = a.cigartuples
                    curr_loc = a.reference_start
                    aln_starts = []
                    aln_ends = []
                    for tup in ct:
                        if tup[0] == 0:
                            aln_starts.append(curr_loc)
                            aln_ends.append(curr_loc + tup[1])
                        curr_loc += tup[1]
                    # If there is any segment that overlaps the peptide without mismatches, accept it.
                    matching_overlap = False
                    overlap = False
                    for e in zip(aln_starts, aln_ends):
                        ol = get_overlap(r[1], r[2], e[0], e[1])
                        if ol:
                            overlap = True
                            # Check if they have mismatches
                            mm = mismatch_in_spliced(a, r, refFasta, suspected_bamchromformat)
                            if mm == 0:
                                matching_overlap = True
                    if matching_overlap:
                        passed_alns.add(a)
                    elif overlap:
                        n_failed_due_to_mismatch += 1
                    else:
                        not_overlapping_peptide += 1
                # Not spliced alignment - just make sure there is no mismatch in the peptide region
                else:
                    unspliced += 1
                    if aln_has_mismatch_in_interval(refFasta, suspected_bamchromformat, a, r[0], int(r[1]), int(r[2])):
                        debug_file.write('Peptide mismatch ' + a.qname + ' ' + str(r[1]) + ' ' + str(r[2]) + '\n')
                        n_failed_due_to_mismatch += 1
                    else:
                        passed_alns.add(a)
        aln_count[p] = len(passed_alns)
        #print('Passed through ' + str(aln_count[p]) + ' alignments in total')
        #print(str(n_failed_due_to_mismatch) + ' rejected due to mismatch')
        #print('Spliced: ' + str(spliced))
        #print('Of these, ' + str(not_overlapping_peptide) + ' did not overlap the peptide segment')
        #print('Unspliced: ' + str(unspliced))
        debug_file.write('BAM file: ' + bam + ' Peptide locus ' + pcoords[0][0] + '-' + str(pcoords[0][1]) + str(pcoords[0][2]) + '\n')
        debug_file.write('Passed through ' + str(len(passed_alns)) + ' alignments. Found ' +  str(n_failed_due_to_mismatch) + ' alignments with mismatches in peptide locus' + '\n#####################\n')
    aln_table[bam]=aln_count

bam = sorted(list(aln_table.keys()))

# Write output file header.
#sys.stdout.write('sequence\ttxtype\tscore\tpsmcount\t')
sys.stdout.write('sequence\t')
for i in range(0,len(bam)):
        sys.stdout.write(bam[i].split('.')[0].split('/')[-1])
        if (i < (len(bam)-1)): sys.stdout.write('\t')
        else: sys.stdout.write('\n')

# And the counts.
for pep in sorted(aln_table[bam[0]].keys()):
    #sys.stdout.write(pep + '\t' + tsv_info['txtype'][pep] + '\t' + tsv_info['score'][pep] + '\t' + tsv_info['psmcount'][pep] + '\t')
    sys.stdout.write(pep + '\t')
    for i in range(0,len(bam)): 
        sys.stdout.write(str(aln_table[bam[i]][pep]))
        if (i < (len(bam)-1)): sys.stdout.write('\t')
    sys.stdout.write('\n')



