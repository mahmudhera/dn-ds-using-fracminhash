from dnds import determine_approx_dnds, determine_correct_dnds

if __name__ == "__main__":
    s1 = 'ACTTGACGAAGCATGCAAGCAACTTGACGAAGCATGCAAGCA'
    s2 = 'ACTTGACGAAGCATGCAAGCATCTTGACGAAGTATGCAAGCA'
    k = 3
    scale_factor = 0.99
    print( determine_approx_dnds(s1, s2, k, scale_factor), determine_correct_dnds(s1, s2) )
