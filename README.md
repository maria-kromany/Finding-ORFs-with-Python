
This project identifies Open Reading Frames (ORFs) in DNA sequences, which are potential gene regions that can code for proteins. ORFs start with "ATG" (start codon) and end with a stop codon (such as "TAA," "TAG," or "TGA"). By locating these segments, we can study gene functions and better understand DNA sequences.

The key tasks include:
Parsing DNA Sequences from FASTA files, using the read_one_seq_fasta function to extract sequence data.
Calculating GC Content via the gc_content function, which computes the percentage of guanine (G) and cytosine (C) bases. This helps filter biologically relevant ORFs.
Identifying ORFs in both DNA strands. Functions like get_orf, one_frame, forward_frames, and reverse_complement work together to locate ORFs across multiple reading frames.
Applying Thresholds with the gene_finder function to filter ORFs based on minimum length and GC content, narrowing down significant gene candidates.


The identified ORFs can be cross-checked with databases like GenBank to uncover biological functions, gene relevance, and human genome presence, connecting computational analysis to real genetic insights.
