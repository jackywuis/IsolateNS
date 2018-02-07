### This program is written for dN/dS ratio calculation based on VCF and GenBank files.

from Bio import SeqIO
from Bio.Seq import Seq
from math import log
from fractions import Fraction
import argparse
import logging

codons = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
          "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
          "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
          "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
          "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
          "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
          "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
          "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
BASES = {'A', 'G', 'T', 'C'}

# copy from https://github.com/adelq/dnds/commit/783e8197541aeb456e19c8ae9effa54cc16b02c0
# copy from https://github.com/adelq/dnds/blob/master/dnds.py
def split_seq(seq, n=3):
    # Returns sequence split into chunks of n characters, default is codons
    return [seq[i:i+n] for i in range(0, len(seq), n)]


def average_list(l1, l2):
    return [(i1 + i2) / 2 for i1, i2 in zip(l1, l2)]


def dna_to_protein(codon):
    # Returns single letter amino acid code for given codon
    return codons[codon]


def translate(seq):
    # Translate a DNA sequence into the 1-letter amino acid sequence
    return "".join([dna_to_protein(codon) for codon in split_seq(seq)])


def is_synonymous(codon1, codon2):
    # Returns boolean whether given codons are synonymous
    return dna_to_protein(codon1) == dna_to_protein(codon2)


def dnds_codon(codon):
    # Returns list of synonymous counts
    syn_list = []
    for i in range(len(codon)):
        base = codon[i]
        other_bases = BASES - {base}
        syn = 0
        for new_base in other_bases:
            new_codon = codon[:i] + new_base + codon[i + 1:]
            syn += int(is_synonymous(codon, new_codon))
        syn_list.append(Fraction(syn, 3))
    return syn_list


def dnds_codon_pair(codon1, codon2):
    # Get the dN/dS for the given codon pair
    return average_list(dnds_codon(codon1), dnds_codon(codon2))


def syn_sum(seq1, seq2):
    # Get the sum of synonymous sites from two DNA sequences
    syn = 0
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        dnds_list = dnds_codon_pair(codon1, codon2)
        syn += sum(dnds_list)
    return syn


def hamming(s1, s2):
    # Return the hamming distance between 2 DNA sequences
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2)) + abs(len(s1) - len(s2))


def codon_subs(codon1, codon2):
    # Returns number of synonymous substitutions in provided codon pair Methodology for multiple substitutions from Dr. Swanson, UWashington https://faculty.washington.edu/wjs18/dnds.ppt
    diff = hamming(codon1, codon2)
    if diff < 1:
        return 0
    elif diff == 1:
        return int(translate(codon1) == translate(codon2))

    syn = 0
    for i in range(len(codon1)):
        base1 = codon1[i]
        base2 = codon2[i]
        if base1 != base2:
            new_codon = codon1[:i] + base2 + codon1[i + 1:]
            syn += int(is_synonymous(codon1, new_codon))
            syn += int(is_synonymous(codon2, new_codon))
    return syn / diff


def substitutions(seq1, seq2):
    # Returns number of synonymous and nonsynonymous substitutions
    dna_changes = hamming(seq1, seq2)
    codon_list1 = split_seq(seq1)
    codon_list2 = split_seq(seq2)
    syn = 0
    for i in range(len(codon_list1)):
        codon1 = codon_list1[i]
        codon2 = codon_list2[i]
        syn += codon_subs(codon1, codon2)
    return (syn, dna_changes - syn)


def clean_sequence(seq):
    # Clean up provided sequence by removing whitespace.
    return seq.replace(' ', '')


def dnds(seq1, seq2):
    """Main function to calculate dN/dS between two DNA sequences per Nei &
    Gojobori 1986. This includes the per site conversion adapted from Jukes &
    Cantor 1967.
    """
    # Strip any whitespace from both strings
    seq1 = clean_sequence(seq1)
    seq2 = clean_sequence(seq2)
    # Check that both sequences have the same length
    assert len(seq1) == len(seq2)
    # Check that sequences are codons
    assert len(seq1) % 3 == 0
    syn_sites = syn_sum(seq1, seq2)
    non_sites = len(seq1) - syn_sites
    logging.info('Sites (syn/nonsyn): {}, {}'.format(syn_sites, non_sites))
    syn_subs, non_subs = substitutions(seq1, seq2)
    pn = non_subs / non_sites
    ps = syn_subs / syn_sites
    dn = -0.75 * log(1 - (4 * pn / 3))
    ds = -0.75 * log(1 - (4 * ps / 3))
    logging.info('dN: {}\t\tdS: {}'.format(round(dn, 3), round(ds, 3)))
    return round(float(pn), 3), round(float(ps), 3), round(float(dn), 3), round(float(ds), 3)


def snp_in_gene(file, record, pos):
    # Confirm whether an SNP in a gene
    output = None
    for records in SeqIO.parse(file, "genbank"):
        if record == records.id:
            for feature in records.features:
                if feature.type == 'CDS' and (feature.location.start <= pos and feature.location.end >= pos):
                    output = feature
                    break
    return output


def get_new_sequence(file, record, input, pos):
    # Predict new sequence where SNP occurs
    for records in SeqIO.parse(file, "genbank"):
        if record == records.id:
            for feature in records.features:
                if 'locus_tag' in feature.qualifiers:
                    if input == feature.qualifiers['locus_tag'][0]:
                        raw_seq = str(feature.extract(records.seq))
                        if feature.location.strand < 0:
                            raw_seq = str(Seq(raw_seq).reverse_complement())
                        new_seq = raw_seq
                        for snp in pos:
                            new_seq = new_seq[: (int(snp[0]) - feature.location.start - 1)] + snp[1] + new_seq[(int(snp[
                                    0]) - feature.location.start):]
                        return raw_seq, new_seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='WooIsolateNS - A tool to calculate . '
                    'Usage: python WooIsolateNS.py -q file.vcf -r reference.gbf '
                    'Output: File_name chromosome pN pS pN/pS dN dS dN/dS'
                    'Copyleft: Jacky Woo from ZHK Research team, iSynBio, SIAT.')
    parser.add_argument('-q', '--query', type=str, help='Your SNP file (in VCF format)')
    parser.add_argument('-r', '--reference', type=str, help='Your reference file (in gbk or gbf format)')
    parser.add_argument('-g', '--gd', action='store_false', required=False,
                        help='Specify that the query file in in .gd format (optional)')
    args = parser.parse_args()
    output = {}
    infh = open(args.query, 'r')
    for line in infh:
        if args.gd:
            if '#' not in line:
                if (len(line.split('\t')[3]) == len(line.split('\t')[4])) and (len(line.split('\t')[4]) == 1):
                    if line.split('\t')[0] not in output:
                        output.update({line.split('\t')[0]: {}})
                    out = snp_in_gene(args.reference, line.split('\t')[0], int(line.split('\t')[1]))
                    if out is not None:
                        out_name = out.qualifiers['locus_tag'][0]
                        if out_name not in output[line.split('\t')[0]]:
                            output[line.split('\t')[0]].update({out_name: [out, []]})
                        output[line.split('\t')[0]][out_name][1] += [
                            (int(line.split('\t')[1]), line.split('\t')[4].split('\n')[0])]
                    out = None
        else:
            if line.split('\t')[0] == 'SNP':
                if line.split('\t')[3] not in output:
                    output.update({line.split('\t')[3]: {}})
                out = snp_in_gene(args.reference, line.split('\t')[3], int(line.split('\t')[4]))
                if out is not None:
                    out_name = out.qualifiers['locus_tag'][0]
                    if out_name not in output[line.split('\t')[3]]:
                        output[line.split('\t')[3]].update({out_name: [out, []]})
                    output[line.split('\t')[3]][out_name][1] += [
                        (int(line.split('\t')[4]), line.split('\t')[5].split('\n')[0])]
                out = None
    infh.close()
    for chromosome in output:
        pn, ps, pns, dn, ds, dns = 0, 0, 0, 0, 0, 0
        for item in sorted(output[chromosome]):
            raw_seq, new_seq = get_new_sequence(args.reference, chromosome, output[chromosome][item][0].qualifiers[
                'locus_tag'][0], output[chromosome][item][1])
            newpn, newps, newdn, newds = dnds(raw_seq, new_seq)
            pn += newpn
            ps += newps
            dn += newdn
            ds += newds
        if pn == 0 and ps == 0:
            pns, dns = '/', '/'
        elif ps == 0:
            pns, dns = '+', '+'
        elif pn == 0:
            pns, dns = '-', '-'
        else:
            pns, dns = pn / ps, dn / ds
        print args.query, chromosome, pn, ps, pns, dn, ds, dns
