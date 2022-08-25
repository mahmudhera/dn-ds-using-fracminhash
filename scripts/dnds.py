from helperfuncs import *

# given two related dna sequences, determine dn/ds
def determine_approx_dnds(dna_seq1, dna_seq2, k, scale_factor):
    nt_containment = sequences_to_containment_using_fmh(dna_seq1, dna_seq2, k, scale_factor)
    aa_containemnt = sequences_to_containment_using_fmh( translate_dna_to_aa(dna_seq1), translate_dna_to_aa(dna_seq2), k, scale_factor)
    dnds = ( 1.0 - aa_containemnt ** (1.0/k) ) / ( aa_containemnt**(1.0/k) - nt_containment**(3.0/k) )
    return dnds
