from helperfuncs import *

# given two related dna sequences, determine dn/ds
def determine_approx_dnds(dna_seq1, dna_seq2, k, scale_factor):
    nt_containment = sequences_to_containment_using_fmh(dna_seq1, dna_seq2, 3*k, scale_factor)
    aa_containemnt = sequences_to_containment_using_fmh( translate_dna_to_aa(dna_seq1), translate_dna_to_aa(dna_seq2), k, scale_factor)
    print(nt_containment, aa_containemnt)
    dnds = ( 1.0 - aa_containemnt ** (1.0/k) ) / ( aa_containemnt**(1.0/k) - nt_containment**(3.0/k) )
    return dnds

def determine_correct_dnds(dna_seq1, dna_seq2):
    num_tota_mutations = 0
    for i in range(len(dna_seq1)):
        if dna_seq1[i] != dna_seq2[i]:
            num_tota_mutations += 1
    aa_seq1 = translate_dna_to_aa(dna_seq1)
    aa_seq2 = translate_dna_to_aa(dna_seq2)
    num_non_syn_mutations = 0
    for i in range(len(aa_seq1)):
        if aa_seq1[i] != aa_seq2[i]:
            num_non_syn_mutations += 1
    print(num_tota_mutations, num_non_syn_mutations)
    return 1.0 * (num_tota_mutations-num_non_syn_mutations)/num_non_syn_mutations
