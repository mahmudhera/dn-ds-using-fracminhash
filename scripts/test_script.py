from dnds import *
from helperfuncs import *
import string
import random
from matplotlib import pyplot as plt

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


    # generate 10k long genome string
    num_runs = 100
    str_len = 10000
    mut_rate = 0.01
    seed = 0
    k = 7
    scale_factor = 0.1
    random.seed(seed)

    correct_dnds_values = []
    estimated_dnds_values_whole_seq = []
    estimated_dnds_values_fmh = []
    for i in range(num_runs):
        s1 = ''.join(random.choices(['A', 'C', 'G', 'T'], k=str_len))
        s2 = mutate_string(s1, mut_rate)

        approx_dnds_using_whole_seq = determine_approx_dnds_no_fmh(s1, s2, k)
        correct_dnds = determine_correct_dnds(s1, s2)
        approx_dnds_using_fmh = determine_approx_dnds_using_fmh(s1, s2, k, scale_factor)

        correct_dnds_values.append(correct_dnds)
        estimated_dnds_values_whole_seq.append(approx_dnds_using_whole_seq)
        estimated_dnds_values_fmh.append(approx_dnds_using_fmh)

    plt.scatter(correct_dnds_values, estimated_dnds_values_whole_seq, label='Using whole seq')
    plt.scatter(correct_dnds_values, estimated_dnds_values_fmh, label='Using FMH with scalef = ' + str(scale_factor))
    tmp = [min(correct_dnds_values), max(correct_dnds_values)]
    plt.plot(tmp, tmp, linestyle='--')
    #plt.ylim(0,100)
    plt.xlabel('Correct dnds value')
    plt.ylabel('Estimated dnds value')
    plt.legend()
    plt.savefig('approx_vs_real_dnds_test.pdf')
