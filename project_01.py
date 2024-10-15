"""
Project 01 Template File

CS/BIOS 112 - Project 01

   Describe
   Project 01 
   Here
   

@author:    <Maria Karim>
UIC NetID:  <mkari5>
Due Date:   <10/10/23>
"""


# The following function should NOT be modified!!
def read_one_seq_fasta(fasta_file):
    """Read a FASTA file that contains one sequence."""
    seq = ''
    with open(fasta_file, 'r') as f:
        f.readline()  # by pass description in first line of data file
        for line in f.readlines():
            seq = seq + line[:-1]
    return seq
# The above function should NOT be modified!!!

def gc_content(seq):
    ''' describe function here '''
    c_count = seq.count("C")
    g_count = seq.count("G")
    total_g_and_c = c_count + g_count
    seq_len = 0 
    for base in seq: 
        seq_len += 1 
    gc_total = total_g_and_c / seq_len 
    return gc_total


# Tests for gc_content function. Should print True in all cases.
print('\ngc_content Tests')
print(gc_content('ATGTGAA') == 0.2857142857142857)
print(gc_content('ATGAGATAAG') == 0.3)


def get_orf(seq):
    ''' describe function here '''
    stop_codons = ['TAG', 'TAA', 'TGA']
    if seq [0:3] == 'ATG' :
        i = 0
        while i < len(seq): 
            cod= seq[i: i + 3]
            if cod in stop_codons: 
                return (seq[0:i])
            i += 3
        return seq
    else: 
        return ''


# Tests for get_orf function. Should print True in all cases.
print('\nget_orf Tests')
print(get_orf('ATGTGAA') == 'ATG')
print(get_orf('ATGAGATAAG') == 'ATGAGA')
print(get_orf('ATGAGATAGG') == 'ATGAGA')
print(get_orf('ATGAGATGAGGGTAA') == 'ATGAGA')
print(get_orf('ATGAAATT') == 'ATGAAATT')


def one_frame(seq):
    ''' describe function here '''
    ocodon = 0
    passlist_codons = []
    x = 0
    while x < len(seq):
        if seq[x:x+3] == "ATG":
            orf = get_orf(seq[x:])
            passlist_codons.append(orf)
            x += len(orf)
        else:
            x += 3
    return passlist_codons

# Tests for one_frame function. Should print True in all cases.
print('\none_frame')
print(one_frame('ATGTGAA') == ['ATG'])
print(one_frame('ATGAGATAAG') == ['ATGAGA'])
print(one_frame('ATGAGATAGG') == ['ATGAGA'])
print(one_frame('ATGAGATGAGGGTAA') == ['ATGAGA'])
print(one_frame('ATGAAATT') == ['ATGAAATT'])
print(one_frame('ATGAGATGAACCATGGGGTAA') == ['ATGAGA', 'ATGGGG'])



def forward_frames(seq):
    ''' describe function here '''
    list_a = one_frame(seq)
    list_b = one_frame(seq[1:])
    list_c = one_frame(seq[2:])
    finalist = list_a + list_b + list_c
    return finalist


# Tests for forward_frames function. Should print True in all cases.
print('\nforward_frames')
print(forward_frames('ATGAGATAAG') == ['ATGAGA'])
print(forward_frames('ATGAGATGAGGGTAA') == ['ATGAGA', 'ATGAGGGTAA'])
print(forward_frames('ATGAAATT') == ['ATGAAATT'])
print(forward_frames('ATGAGATGACACCATGGGGTAA') == ['ATGAGA', 'ATGGGG', 'ATGACACCATGGGGTAA'])


def reverse_complement(seq):
    '''describe function here ''' 

    reverse_complement = seq[::-1]
    x = "" 
    i = 0
    for i in range(len(reverse_complement)):
        if reverse_complement[i] == "A":
            x = x + "T"
        elif reverse_complement[i] == 'T':
            x = x + "A" 
        elif reverse_complement[i] == "G": 
            x = x + "C" 
        elif reverse_complement[i] == "C": 
            x = x + "G"
    return x
# Tests for reverse_complement function. Should print True.
print ("\nreverse_compliment")
print (reverse_complement('ATGCTTG') == 'CAAGCAT')
print (reverse_complement('AAAGGG') == 'CCCTTT')
print (reverse_complement('TTTCCC') == 'GGGAAA')
print (reverse_complement('ATCGATCAGTCCTAGCATCG') == 'CGATGCTAGGACTGATCGAT')

def gene_finder(fasta_file, min_len, min_gc):
    ''' describe function here '''
    find_dna = read_one_seq_fasta(fasta_file)
    lista = forward_frames(find_dna)
    revseq = reverse_complement(find_dna)
    listb = forward_frames(revseq) 
    listc = lista + listb 
    i = 0 
    final = [] 
    
    for orf in listc: 
        if len(orf) >= min_len and gc_content(orf) >= min_gc: 
            final.append([len(orf), gc_content(orf), orf]) 
    return final


# Tests for gene_finder function. Should print True.
print('\ngene_finder')
calculated_result = gene_finder('gene_finder_test.fasta', 6, 0.45)
desired_result = [[6, 0.6666666666666666, 'ATGCCC'], [9, 0.7777777777777778, 'ATGCCCCGG']]
print( calculated_result == desired_result)

orf_list = gene_finder('human_chr9_segment.fasta', 550, 0.45)
print (len(orf_list) == 5)

# viewing the results of the gene_finder( ) calculations
orf_list = gene_finder('gene_finder_test.fasta', 6, 0.45)
print (orf_list)


"""
Identify ORFs in either of the provided fasta files

< In this comment, replace this paragraph with the information you learn from the 
  top BLAST hits when using data produced by the results of gene_finder( ) when using 
  the data in the X73525.fasta or the human chromosome X partial fasta file on the GenBank website at:
       http://www.ncbi.nlm.nih.gov/blast/Blast.cgi
  Include the parameter values used when you called gene_finder( ).  
  Briefly describe the likely function of your gene and paste your gene information 
  and sequence (length, %GC, and DNA sequence) that you used in the BLAST search.

  Be sure sure to include the following information:
  a. Why is this gene particularly relevant today? (use gene from file: human chromosome X partial.fasta)
  b. Do YOU have this gene? (use gene from file: human chromosome X partial.fasta)
  c. How many ORFs you identified in the sequence? (indicate which fasta file)
  d. What is the ORFs function? (indicate which fasta file)
  >
"""