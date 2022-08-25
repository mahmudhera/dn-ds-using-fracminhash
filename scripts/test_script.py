from dnds import determine_approx_dnds

if __name__ == "__main__":
    s1 = 'ACTTGACGAAGCATGCAAGCA'
    s2 = 'ACTAGAAGAAGCATACTAGCA'
    k = 7
    scale_factor = 1.0
    print( determine_approx_dnds(s1, s2, k, scale_factor) )
