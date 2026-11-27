#!/usr/bin/env python3
#
# MEIsensor annotation, computed by MEIsensor itself from the reference genome:
#   - EN motif:      6-bp T-tract/A-tract core at the physical TSD boundary (-6..+6 window)
#   - TSD:           exact (0-mismatch) target-site duplication, internal implementation
#   - poly(A/T):     terminal A/T-rich tract of the inserted sequence after removing the TSD
#
# No external annotation tool is called and no third-party annotation is imported. The
# TSD/poly(A)/EN fields are written into the output VCF INFO column by vcf.py (--annotate).
#

import logging

log = logging.getLogger('MEIsensor.annotate')

#
# EN motif annotation
#
# The EN nick is anchored at the physical TSD boundary, not found by scanning the
# reference for a T/A junction:
#   - POLYA_STRAND = "+": first nick on the bottom strand at the TSD 5' boundary;
#     the 12-bp window (positions -6..-1, +1..+6) is rc(ref[start-6 : start+6]).
#   - POLYA_STRAND = "-": first nick on the top strand at the TSD 3' boundary;
#     the window is ref[end-6 : end+6].
# EN_MOTIF reports the 6-bp core (positions -4..-1 and +1..+2) as a T-tract/A-tract
# pair; EN_STRAND is the nicked strand and EN_OFFSET is 0 because the junction is
# defined by the TSD boundary itself.
#
# The TSD interval is derived in annotate_insertion (prefix: start = POS,
# end = POS + TSD_LEN; suffix: start = POS - TSD_LEN, end = POS).

EN_FLANK = 6          # 6 bp on each side of the nick (12-bp window)
EN_CORE_START = 2     # window index of position -4
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


#
# Target-site duplication (TSD)
#
# Definition: the TSD is the sequence of the target site that is duplicated on both
# sides of the insertion. It is measured as the longest EXACT (0-mismatch) match
# between the inserted sequence and the reference sequence immediately flanking
# the insertion point:
#
#   (A) breakpoint before the target site: insert starts with the TSD copy
#       -> compare the insert prefix with the reference bases AFTER POS;
#   (B) breakpoint after the target site:  insert ends with the TSD copy
#       -> compare the insert suffix with the reference bases ending AT POS.
#
# Only one comparison is used (the longer of the two; they are never summed), no
# mismatches are tolerated and the search stops at the first difference. A TSD is
# called only when the exact match is at least MIN_TSD bp; otherwise the length is
# reported as 0. The VCF anchor base is not part of the inserted sequence and is
# not counted.

MIN_TSD = 4        # minimum exact match length reported as a TSD (bp)
TSD_MAX_LEN = 60   # compare at most this many bp


def exact_prefix(a, b, max_len=TSD_MAX_LEN):
    """
    Length of the exact (0-mismatch) common prefix of a and b, capped at max_len.
    Stops at the first mismatching base.
    """
    n = 0
    for x, y in zip(a[:max_len], b[:max_len]):
        if x != y:
            break
        n += 1
    return n


#
# poly(A/T) tail
#
# Parameters follow the conventional 10-bp-window monomer rule used for poly(A/T)
# calling in MEI annotation (10-bp windows, >=80% purity, monomer >=10 bp), so that
# the lengths are directly comparable with other tools' poly(A) fields. The
# sequence handed to the search is the inserted sequence with the exact TSD
# removed. A tract is only accepted when every consecutive window at the terminus
# passes the purity threshold (POLYA_MAX_WIN_DIST = 0: gaps of non-monomeric
# windows are not bridged); this is the definition used for the reported
# distributions.

POLYA_WINDOW = 10        # bp
POLYA_MIN_MONOMER = 10   # bp
POLYA_MIN_PURITY = 80    # percent of the target base per window
POLYA_MAX_WIN_DIST = 0   # number of non-monomeric windows tolerated (0 = none)


def find_poly_tract(seq, base, from_3prime=True):
    """
    Length of the terminal tract of `base` in seq.

    Input:
        1. seq: sequence to scan (upper case)
        2. base: 'A' (poly(A) at the 3' end) or 'T' (poly(T) at the 5' end)
        3. from_3prime: scan from the 3' end (True) or from the 5' end (False)

    Output:
        tract length in bp (0 when no tract of at least POLYA_MIN_MONOMER is found)
    """
    n = len(seq)
    if n < POLYA_MIN_MONOMER:
        return 0

    need = POLYA_MIN_PURITY * POLYA_WINDOW / 100.0
    length = 0
    if from_3prime:
        for start in range(n - POLYA_WINDOW, -1, -1):
            window = seq[start:start + POLYA_WINDOW]
            if window.count(base) >= need:
                length = n - start
            else:
                break
    else:
        for end in range(POLYA_WINDOW, n + 1):
            window = seq[end - POLYA_WINDOW:end]
            if window.count(base) >= need:
                length = end
            else:
                break
    return length if length >= POLYA_MIN_MONOMER else 0


def annotate_insertion(reference_handle, contig, pos, inserted_seq):
    """
    Compute the strict TSD, poly(A/T) and TSD-anchored EN annotation of one insertion.

    Input:
        1. reference_handle: open pysam.FastaFile
        2. contig, pos: VCF CHROM and 1-based POS of the base preceding the insertion
        3. inserted_seq: the inserted sequence alone (the VCF anchor base is NOT part of it;
           the written ALT is anchor + inserted_seq)

    Output:
        dict with keys tsd_len, tsd_seq, polya_len, polya_seq, polya_strand,
        en_motif, en_strand, en_offset, or None when the call cannot be annotated
        (empty ALT, missing contig). EN fields are empty/0 when no exact TSD
        (>= MIN_TSD) or no poly(A) strand is available, because the EN window is
        anchored at the physical TSD boundary.
    """
    if reference_handle is None or not inserted_seq or inserted_seq.startswith('<'):
        return None

    insert = inserted_seq.upper()
    if not insert:
        return None

    try:
        ref_right = reference_handle.fetch(contig, pos, pos + TSD_MAX_LEN).upper()
        ref_left = reference_handle.fetch(contig, max(0, pos - TSD_MAX_LEN), pos).upper()
    except (KeyError, ValueError):
        return None

    # (A) insert prefix vs reference after POS ; (B) insert suffix vs reference ending at POS
    len_a = exact_prefix(insert, ref_right)
    len_b = exact_prefix(insert[::-1], ref_left[::-1])
    if len_a >= len_b:
        tsd_len, tsd_seq, core = len_a, insert[:len_a], insert[len_a:]
        side = 'right'
    else:
        tsd_len, tsd_seq, core = len_b, insert[-len_b:] if len_b else '', insert[:-len_b] if len_b else insert
        side = 'left'

    if tsd_len < MIN_TSD:
        tsd_len, tsd_seq, core = 0, '', insert

    polya_len, polya_seq, polya_strand = 0, '', ''
    if core:
        a_len = find_poly_tract(core, 'A', from_3prime=True)
        t_len = find_poly_tract(core, 'T', from_3prime=False)
        if a_len or t_len:
            if a_len >= t_len:
                polya_len, polya_seq, polya_strand = a_len, core[-a_len:], '+'
            else:
                polya_len, polya_seq, polya_strand = t_len, core[:t_len], '-'

    # EN motif: fixed 12-bp window (-6..-1, +1..+6) anchored at the physical TSD
    # boundary; the 6-bp core (-4..-1, +1, +2) is reported as T-tract/A-tract.
    en_motif, en_strand, en_offset = '', '', 0
    if tsd_len >= MIN_TSD and polya_strand in ('+', '-'):
        start = pos if side == 'right' else pos - tsd_len
        end = pos + tsd_len if side == 'right' else pos
        try:
            if polya_strand == '+':
                en_window = reference_handle.fetch(contig, start - EN_FLANK, start + EN_FLANK).upper()
                en_window = en_window.translate(_COMPLEMENT)[::-1]
                en_strand = '-'
            else:
                en_window = reference_handle.fetch(contig, end - EN_FLANK, end + EN_FLANK).upper()
                en_strand = '+'
        except (KeyError, ValueError):
            en_window = ''
        if len(en_window) == 2 * EN_FLANK and not set(en_window) - set('ACGT'):
            core6 = en_window[EN_CORE_START:EN_CORE_START + 6]
            en_motif = core6[:4] + '/' + core6[4:]
            en_offset = 0
        else:
            en_motif, en_strand = '', ''

    return {
        'tsd_len': tsd_len,
        'tsd_seq': tsd_seq,
        'polya_len': polya_len,
        'polya_seq': polya_seq,
        'polya_strand': polya_strand,
        'en_motif': en_motif,
        'en_strand': en_strand,
        'en_offset': en_offset,
    }
