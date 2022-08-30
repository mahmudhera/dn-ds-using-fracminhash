from dnds import *
from helperfuncs import *

if __name__ == "__main__":
    s1 = 'ACTTGACGAAGGATGCAAGCAACTTGACGAAGCATGCAAGCT'
    s2 = 'ACCTGACGAAGCATGCAAGCAACTTGACGAAGCATGCAAGCA'
    k = 7
    scale_factor = 1.0

    tralated_s1 = translate_dna_to_aa(s1)
    print(len(tralated_s1), tralated_s1)

    T = translate_dna_to_aa
    aa_mut_rate = sequences_to_mut_rate_perfect(T(s1), T(s2), 3)
    print(aa_mut_rate)

    approx_dnds = determine_approx_dnds_no_fmh(s1, s2, k)
    print(approx_dnds)

    correct_dnds = determine_correct_dnds(s1, s2)
    print(correct_dnds)
