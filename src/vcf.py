#!/usr/bin/env python3
#
# MEIsensor

import logging

import pysam
import os

import SV
import utils
import annotate
from config import Config
from trained_model.model import classify_INS


log = logging.getLogger(__name__)


def format_info(k, v):
    if isinstance(v, float):
        return f"{k}={v:.3f}"
    elif isinstance(v, list):
        return f"{k}={','.join(v)}"
    else:
        return f"{k}={v}"




def format_genotype(gt):
    """
    Format a genotype tuple as VCF sample column value.
    Accepts a 5-tuple (a, b, qual, dr, dv) or a 6-tuple (a, b, qual, dr, dv, svid).
    """
    if len(gt) == 5:
        a, b, qual, dr, dv = gt
        gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}"
    else:
        a, b, qual, dr, dv, svid = gt
        gt_sep = "/"
        return f"{a}{gt_sep}{b}:{qual}:{dr}:{dv}:{svid}"



class VCF:
    def __init__(self, config: Config, handle):
        self.config = config
        self.handle = handle
        self.call_count = 0
        self.info_order = ["SVTYPE", "SVLEN", "END", "SUPPORT", "RNAMES", "COVERAGE", "STRAND"]
        if config.qc_nm_measure:
            self.info_order.append("NM")

        if config.dev_emit_sv_lengths:
            self.info_order.append("SVLENGTHS")

        self.default_genotype = config.genotype_none

        # Add phasing if needed
        self.genotype_format = config.genotype_format
        self.reference_handle = None
        self.header_str = ""

    def open_reference(self):
        if self.config.reference is None:
            return

        if not os.path.exists(self.config.reference + ".fai") and not os.path.exists(self.config.reference + ".gzi"):
            print(f"[MEIsensor] No index found for {self.config.reference}; building it with pysam.faidx "
                  f"(this can take a while)")
            pysam.faidx(self.config.reference)
        self.reference_handle = pysam.FastaFile(self.config.reference)

    def write_header(self, contigs_lengths):
        self.write_header_line("fileformat=VCFv4.2")
        self.write_header_line(f"source=MEIsensor-{self.config.build}")
        for contig, contig_len in contigs_lengths:
            self.write_header_line(f"contig=<ID={contig},length={contig_len}>")

        self.write_header_line('ALT=<ID=INS,Description="Mobile element insertion">')
        self.write_header_line('FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        self.write_header_line('FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        self.write_header_line('FORMAT=<ID=DR,Number=1,Type=Integer,Description="Reads supporting the reference allele">')
        self.write_header_line('FORMAT=<ID=DV,Number=1,Type=Integer,Description="Reads supporting the insertion">')

        self.write_header_line('FILTER=<ID=PASS,Description="Passed all filters">')
        self.write_header_line('FILTER=<ID=GT,Description="Genotype confidence below threshold">')
        self.write_header_line('FILTER=<ID=SUPPORT_MIN,Description="Supporting read count below minimum">')
        self.write_header_line('FILTER=<ID=STDEV_POS,Description="Insertion point spread too large">')
        self.write_header_line('FILTER=<ID=STDEV_LEN,Description="Insertion length spread too large">')
        self.write_header_line('FILTER=<ID=COV_MIN_GT,Description="Coverage too low for reliable genotyping">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_INS,Description="Coverage drops around the insertion">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_US,Description="Coverage change between upstream and start">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_SC,Description="Coverage change between start and center">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_CE,Description="Coverage change between center and end">')
        self.write_header_line('FILTER=<ID=COV_CHANGE_FRAC_ED,Description="Coverage change between end and downstream">')
        self.write_header_line('FILTER=<ID=MOSAIC_VAF,Description="Allele fraction outside the mosaic range">')
        self.write_header_line('FILTER=<ID=NOT_MOSAIC_VAF,Description="Allele fraction outside the non-mosaic range">')
        self.write_header_line('FILTER=<ID=ALN_NM,Description="Too many mismatches in supporting alignments">')
        self.write_header_line('FILTER=<ID=STRAND,Description="Supported by one strand only">')
        self.write_header_line('FILTER=<ID=STRAND_MOSAIC,Description="Mosaic call supported by one strand only">')
        self.write_header_line('FILTER=<ID=SVLEN_MIN,Description="Insertion shorter than the minimum length">')
        self.write_header_line('INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Insertion point resolved within the precision window">')
        self.write_header_line('INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Insertion point imprecise">')
        self.write_header_line('INFO=<ID=MOSAIC,Number=0,Type=Flag,Description="Putative mosaic insertion">')
        self.write_header_line('INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the inserted sequence">')
        self.write_header_line('INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the insertion event">')
        self.write_header_line('INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="Number of reads supporting the insertion">')
        self.write_header_line('INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="Soft-clipped reads suggesting a long insertion">')
        self.write_header_line('INFO=<ID=END,Number=1,Type=Integer,Description="End position of the event on the reference">')
        self.write_header_line('INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="Spread of insertion points across supporting reads">')
        self.write_header_line('INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="Spread of insertion lengths across supporting reads">')
        self.write_header_line('INFO=<ID=COVERAGE,Number=.,Type=Float,Description="Depth around the insertion: upstream, start, center, end, downstream">')
        self.write_header_line('INFO=<ID=STRAND,Number=1,Type=String,Description="Strands of the supporting reads">')
        self.write_header_line('INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count summed over samples">')
        self.write_header_line('INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Per-sample support vector">')
        self.write_header_line('INFO=<ID=RNAMES,Number=.,Type=String,Description="Names of supporting reads, when enabled">')
        self.write_header_line('INFO=<ID=VAF,Number=1,Type=Float,Description="Fraction of supporting reads">')
        self.write_header_line('INFO=<ID=NM,Number=.,Type=Float,Description="Mismatch rate of supporting alignments">')
        if self.config.annotate:
            self.write_header_line('INFO=<ID=EN_MOTIF,Number=1,Type=String,Description="EN cleavage motif at the '
                                   'TSD boundary: 6-bp core (positions -4..-1 and +1..+2) reported as '
                                   'T-tract/A-tract on the nicked strand; absent when no exact TSD (>=4 bp) or '
                                   'no poly(A) strand is available">')
            self.write_header_line('INFO=<ID=EN_STRAND,Number=1,Type=String,Description="Strand carrying the EN '
                                   'nick at the TSD boundary (+ or -)">')
            self.write_header_line('INFO=<ID=EN_OFFSET,Number=1,Type=Integer,Description="Offset in bp between '
                                   'the EN junction and the TSD boundary (0 under the TSD-anchored definition)">')
            self.write_header_line('INFO=<ID=TSD_LEN,Number=1,Type=Integer,Description="Target-site duplication '
                                   'length: longest exact (0-mismatch) match between the inserted sequence and the '
                                   'reference immediately flanking the insertion, requiring at least 4 bp; 0 when '
                                   'no duplication is detected">')
            self.write_header_line('INFO=<ID=TSD_SEQ,Number=1,Type=String,Description="Target-site duplication '
                                   'sequence (exact match, no mismatches); absent when TSD_LEN=0">')
            self.write_header_line('INFO=<ID=POLYA_LEN,Number=1,Type=Integer,Description="poly(A) tail length of '
                                   'the inserted sequence after removal of the TSD: terminal A-rich (3\' end) or '
                                   'T-rich (5\' end) tract of at least 10 bp with at least 80% purity in every '
                                   '10-bp window; 0 when no tail is detected">')
            self.write_header_line('INFO=<ID=POLYA_SEQ,Number=1,Type=String,Description="poly(A/T) tract sequence '
                                   'as found in the inserted sequence; absent when POLYA_LEN=0">')
            self.write_header_line('INFO=<ID=POLYA_STRAND,Number=1,Type=String,Description="Orientation of the '
                                   'poly(A) tract: + if A-rich at the 3\' end of the inserted sequence, - if '
                                   'T-rich at the 5\' end">')
        samples_header = "\t".join(sample_id for _, sample_id in self.config.sample_ids_vcf)
        self.write_raw(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples_header}")
        
    def write_raw(self, text, endl="\n"):
        if self.config.vcf_output_bgz:
            self.handle.write(text.encode())
            self.handle.write(endl.encode())
        else:
            self.handle.write(text)
            self.handle.write(endl)

    def write_header_line(self, text):
        self.write_raw("##" + text)

    def write_call(self, call):
        if call.svtype != "INS":
            return
        if call.alt == "<INS>":
            return
        # pysam coordinates are 0-based, VCF 1-based
        # but VCF also requires the index of the base before the SV to be reported,
        # so we are fine without offsetting
        end = call.end
        pos = call.pos if call.pos > 0 else 1

        # Determine genotypes columns
        ac = 0
        supvec = []
        sample_genotypes = []
        for internal_id, _ in self.config.sample_ids_vcf:
            if internal_id in call.genotypes and call.genotypes[internal_id] is not None:
                gt_curr = call.genotypes[internal_id]
                sample_genotypes.append(format_genotype(gt_curr))
                if gt_curr[0] != "." and gt_curr[4] > 0:  # Not non-genotype and has supporting reads
                    ac += sum(call.genotypes[internal_id][:2])
                    supp = "1"
                else:
                    supp = "0"
            else:
                sample_genotypes.append(format_genotype(self.default_genotype))
                supp = "0"
            supvec.append(supp)

        if len(self.config.sample_ids_vcf) > 1:
            call.set_info("AC", ac)
            call.set_info("SUPP_VEC", svec := "".join(supvec))

            if int(svec) == 0:
                log.debug(f'Skipped {call}: no sample has support.')
                return

            if ac == 0:
                call.filter = "GT"

        if (self.config.annotate and call.svtype == "INS" and call.alt != '<INS>'
                and self.reference_handle is not None):
            # Strict TSD, poly(A/T) and TSD-anchored EN motif annotation, computed by
            # MEIsensor itself from the reference genome and the inserted sequence
            # (no external annotation tool). call.alt is the inserted sequence; the
            # written ALT is call.ref + call.alt.
            hit = annotate.annotate_insertion(self.reference_handle, call.contig, pos, call.alt)
            if hit is not None:
                call.set_info("TSD_LEN", hit['tsd_len'])
                if hit['tsd_len']:
                    call.set_info("TSD_SEQ", hit['tsd_seq'])
                call.set_info("POLYA_LEN", hit['polya_len'])
                if hit['polya_len']:
                    call.set_info("POLYA_SEQ", hit['polya_seq'])
                    call.set_info("POLYA_STRAND", hit['polya_strand'])
                if hit['en_motif']:
                    call.set_info("EN_MOTIF", hit['en_motif'])
                    call.set_info("EN_STRAND", hit['en_strand'])
                    call.set_info("EN_OFFSET", hit['en_offset'])

        # Output core SV attributes
        infos = {
            "SVTYPE": call.svtype,
            "SVLEN": call.svlen,
            "SVLENGTHS": ",".join(map(str, call.svlens)) if call.svlens else None,
            "END": end,
            "SUPPORT": call.support,
            "RNAMES": call.rnames if self.config.output_rnames else None,
            "COVERAGE": f"{call.coverage_upstream},{call.coverage_start},{call.coverage_center},{call.coverage_end},"
                        f"{call.coverage_downstream}",
            "STRAND": ("+" if call.fwd > 0 else "") + ("-" if call.rev > 0 else ""),
            "NM": call.nm
        }
        infos_ordered = ["PRECISE" if call.precise else "IMPRECISE"]
        af = call.get_info("VAF")
        af = af if af is not None else 0
        sv_is_mosaic = af <= self.config.mosaic_af_max
        if sv_is_mosaic and self.config.mosaic:
            infos_ordered.append("MOSAIC")
        infos_ordered.extend(format_info(k, infos[k]) for k in self.info_order if infos[k] is not None)
        info_str = ";".join(infos_ordered)

        # Output call specific additional information
        for k in sorted(call.info):
            if call.info[k] is None:
                continue
            info_str += ";" + format_info(k, call.info[k])

        # if call.id==None:
        #    call.id=f"{call.svtype}.{self.call_count+1:06}"

        if (not self.config.symbolic and call.svtype == "DEL" and self.reference_handle is not None
                and abs(call.svlen) <= self.config.max_del_seq_len):
            return

        if self.config.symbolic:
            call.ref = "N"
            call.alt = f"<{call.svtype}>"
        else:
            if self.reference_handle is not None and call.ref == 'N':
                try:
                    call.ref = self.reference_handle.fetch(call.contig, start := max(0, call.pos - 1), start + 1)
                except (KeyError, ValueError):
                    ...
                else:
                    if call.svtype == "INS" and call.alt != '<INS>':
                        INS_seq = call.ref + call.alt
                        Alu_INS, LINE1_INS, SVA_INS, other_INS, score = classify_INS(INS_seq, self.config)
                        if SVA_INS or Alu_INS or LINE1_INS:
                            call.alt = INS_seq
                            if SVA_INS:
                                call.id = 'SVA.INS'
                            elif Alu_INS:
                                call.id = "ALU.INS"
                            elif LINE1_INS:
                                call.id = f"LINE1.INS"
                        else:
                            log.debug(f"Skipped {call}: not an Alu/LINE1/SVA insertion.")
                            return

        call.qual = max(0, min(60, call.qual)) if call.qual is not None else None

        self.write_raw("\t".join(str(v) for v in [call.contig, pos, self.config.id_prefix + call.id, call.ref,
                                                  call.alt, call.qual if call.qual is not None else '.', call.filter, info_str, self.genotype_format] +
                                 sample_genotypes))
        self.call_count += 1

    def read_svs_iter(self):
        self.header_str = ""
        line_index = 0
        for line in self.handle:
            try:
                if isinstance(line, bytes):
                    line = line.decode("utf-8")
                line_index += 1
                line_strip = line.strip()
                if line_strip == "" or line_strip[0] == "#":
                    if line_strip[0] == "#":
                        self.header_str += line_strip + "\n"
                    continue
                CHROM, POS, _, REF, ALT, QUAL, FILTER, INFO = line.split("\t")[:8]
                info_dict = {}
                for info_item in INFO.split(";"):
                    if "=" in info_item:
                        key, value = info_item.split("=")
                    else:
                        key, value = info_item, True
                    info_dict[key] = value
                call = SV.SVCall(contig=CHROM,
                                 pos=int(POS) - 1,
                                 id=line_index,
                                 ref=REF,
                                 alt=ALT,
                                 qual=int(QUAL) if QUAL != '.' else None,
                                 filter=FILTER,
                                 info=info_dict,
                                 svtype=None,
                                 svlen=None,
                                 end=None,
                                 rnames=None,
                                 qc=True,
                                 postprocess=None,
                                 genotypes=None,
                                 precise=None,
                                 support=0,
                                 fwd=0,
                                 rev=0,
                                 nm=-1)
                if len(call.alt) > len(call.ref):
                    call.svtype = "INS"
                    call.svlen = len(call.alt)
                    call.end = call.pos
                else:
                    call.svtype = "DEL"
                    call.svlen = -len(call.ref)
                    call.end = call.pos + call.svlen

                if "SVTYPE" in info_dict:
                    call.svtype = info_dict["SVTYPE"]
                    if call.svtype == "TRA":
                        call.svtype = "BND"

                if "SVLEN" in info_dict:
                    call.svlen = int(info_dict["SVLEN"])
                if "SVLENGTHS`" in info_dict:
                    call.svlens = info_dict["SVLENGTHS"]

                if "END" in info_dict:
                    call.end = int(info_dict["END"])

                call.raw_vcf_line = line_strip
                call.raw_vcf_line_index = line_index
                yield call
            except Exception as e:
                utils.fatal_error(f"Error parsing input VCF: Line {line_index}: {e}")

    def rewrite_genotype(self, svcall):
        parts_no_gt = svcall.raw_vcf_line.split("\t")[:8]
        gt_format = self.config.genotype_format
        if svcall.genotype_match_sv != None:
            if len(svcall.genotype_match_sv.genotypes) > 0:
                gt = svcall.genotype_match_sv.genotypes[0]
            else:
                gt = svcall.genotypes[0]
        else:
            gt = svcall.genotypes[0]
        parts = parts_no_gt + [gt_format, format_genotype(gt)]
        self.write_raw("\t".join(parts))

    def rewrite_header_genotype(self,orig_header):
        header_lines=orig_header.split("\n")
        header_lines.insert(1,'##genotypeDate="'+self.config.start_date+'"')
        header_lines.insert(1,'##genotypeCommand="'+self.config.command+'"')
        header_lines.insert(1,f"##genotypeSource=MEIsensor-{self.config.build}")

        has_gt_headers = {
            "GT": False,
            "GQ": False,
            "DR": False,
            "DV": False,
        }
        for header_line in header_lines:
            for gt in has_gt_headers.keys():
                if "##FORMAT=<ID="+gt+"," in header_line:
                    has_gt_headers[gt] = True

        if not has_gt_headers["GT"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        if not has_gt_headers["GQ"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
        if not has_gt_headers["DR"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">')
        if not has_gt_headers["DV"]:
            header_lines.insert(len(header_lines)-2, '##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">')

        self.write_raw("\n".join(header_lines), endl="")

    def close(self):
        self.handle.close()