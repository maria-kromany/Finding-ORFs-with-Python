"""
Project 01 Template File

CS/BIOS 112 - Project 01

This project involves analyzing DNA sequences to identify Open Reading Frames (ORFs)
that could represent potential genes. The code includes functions to calculate GC content,
find ORFs, and retrieve the reverse complement of a DNA sequence, all of which are
key in gene prediction and analysis.

@author:    <Maria Karim>
UIC NetID:  <mkari5>
Due Date:   <10/10/23>
"""

# Function to read a single sequence from a FASTA file
# The following function should NOT be modified!!
def read_one_seq_fasta(fasta_file):
    """Read a FASTA file that contains one sequence."""
    seq = ''
    with open(fasta_file, 'r') as f:
        f.readline()  # Skip the description line in the FASTA file
        for line in f.readlines():
            seq = seq + line[:-1]  # Append each line to the sequence, removing newline characters
    return seq
# The above function should NOT be modified!!!

# Function to calculate GC content in a DNA sequence
def gc_content(seq):
    '''Calculate the GC content of a DNA sequence as a ratio.'''
    c_count = seq.count("C")  # Count cytosine bases in the sequence
    g_count = seq.count("G")  # Count guanine bases in the sequence
    total_g_and_c = c_count + g_count  # Total GC content
    seq_len = 0  # Initialize sequence length
    for base in seq: 
        seq_len += 1  # Calculate total sequence length
    gc_total = total_g_and_c / seq_len  # Calculate GC content as a fraction of total length
    return gc_total

# Tests for gc_content function. Should print True in all cases.
print('\ngc_content Tests')
print(gc_content('ATGTGAA') == 0.2857142857142857)
print(gc_content('ATGAGATAAG') == 0.3)

# Function to find an ORF starting with 'ATG' and ending at a stop codon
def get_orf(seq):
    '''Find the longest ORF in a DNA sequence starting with ATG and ending at a stop codon.'''
    stop_codons = ['TAG', 'TAA', 'TGA']  # List of stop codons
    if seq[0:3] == 'ATG':  # Check if sequence starts with a start codon 'ATG'
        i = 0
        while i < len(seq): 
            cod = seq[i: i + 3]  # Extract each codon from the sequence
            if cod in stop_codons:  # Stop if a stop codon is found
                return seq[0:i]  # Return ORF up to the stop codon
            i += 3
        return seq  # If no stop codon, return entire sequence
    else: 
        return ''  # Return empty if no 'ATG' at start

# Tests for get_orf function. Should print True in all cases.
print('\nget_orf Tests')
print(get_orf('ATGTGAA') == 'ATG')
print(get_orf('ATGAGATAAG') == 'ATGAGA')
print(get_orf('ATGAGATAGG') == 'ATGAGA')
print(get_orf('ATGAGATGAGGGTAA') == 'ATGAGA')
print(get_orf('ATGAAATT') == 'ATGAAATT')

# Function to find all ORFs in one reading frame
def one_frame(seq):
    '''Find all ORFs in a single reading frame of a DNA sequence.'''
    passlist_codons = []
    x = 0
    while x < len(seq):
        if seq[x:x+3] == "ATG":  # Check if start codon 'ATG' is found
            orf = get_orf(seq[x:])  # Get the ORF starting at 'ATG'
            passlist_codons.append(orf)  # Append the ORF to the list
            x += len(orf)  # Move index forward by length of ORF
        else:
            x += 3  # Move index forward by one codon if no 'ATG'
    return passlist_codons

# Tests for one_frame function. Should print True in all cases.
print('\none_frame')
print(one_frame('ATGTGAA') == ['ATG'])
print(one_frame('ATGAGATAAG') == ['ATGAGA'])
print(one_frame('ATGAGATAGG') == ['ATGAGA'])
print(one_frame('ATGAGATGAGGGTAA') == ['ATGAGA'])
print(one_frame('ATGAAATT') == ['ATGAAATT'])
print(one_frame('ATGAGATGAACCATGGGGTAA') == ['ATGAGA', 'ATGGGG'])

# Function to find ORFs in all forward frames (0, 1, and 2) of the sequence
def forward_frames(seq):
    '''Find ORFs in all three forward reading frames of a DNA sequence.'''
    list_a = one_frame(seq)       # ORFs in frame 0
    list_b = one_frame(seq[1:])    # ORFs in frame 1
    list_c = one_frame(seq[2:])    # ORFs in frame 2
    finalist = list_a + list_b + list_c  # Combine ORFs from all frames
    return finalist

# Tests for forward_frames function. Should print True in all cases.
print('\nforward_frames')
print(forward_frames('ATGAGATAAG') == ['ATGAGA'])
print(forward_frames('ATGAGATGAGGGTAA') == ['ATGAGA', 'ATGAGGGTAA'])
print(forward_frames('ATGAAATT') == ['ATGAAATT'])
print(forward_frames('ATGAGATGACACCATGGGGTAA') == ['ATGAGA', 'ATGGGG', 'ATGACACCATGGGGTAA'])

# Function to find the reverse complement of a DNA sequence
def reverse_complement(seq):
    '''Return the reverse complement of a DNA sequence.'''
    reverse_complement = seq[::-1]  # Reverse the DNA sequence
    x = "" 
    for base in reverse_complement:
        if base == "A":
            x = x + "T"  # Replace 'A' with 'T'
        elif base == 'T':
            x = x + "A"  # Replace 'T' with 'A'
        elif base == "G": 
            x = x + "C"  # Replace 'G' with 'C'
        elif base == "C": 
            x = x + "G"  # Replace 'C' with 'G'
    return x

# Tests for reverse_complement function. Should print True.
print ("\nreverse_compliment")
print (reverse_complement('ATGCTTG') == 'CAAGCAT')
print (reverse_complement('AAAGGG') == 'CCCTTT')
print (reverse_complement('TTTCCC') == 'GGGAAA')
print (reverse_complement('ATCGATCAGTCCTAGCATCG') == 'CGATGCTAGGACTGATCGAT')

# Function to find potential genes based on ORF length and GC content criteria
def gene_finder(fasta_file, min_len, min_gc):
    '''Identify ORFs in a DNA sequence that meet minimum length and GC content requirements.'''
    find_dna = read_one_seq_fasta(fasta_file)  # Read DNA sequence from FASTA file
    lista = forward_frames(find_dna)           # Find ORFs in forward frames
    revseq = reverse_complement(find_dna)      # Reverse complement of DNA sequence
    listb = forward_frames(revseq)             # Find ORFs in reverse frames
    listc = lista + listb                      # Combine ORFs from both directions
    final = [] 
    
    # Filter ORFs based on minimum length and GC content
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

# Viewing the results of the gene_finder calculations
orf_list = gene_finder('gene_finder_test.fasta', 6, 0.45)
print (orf_list)
