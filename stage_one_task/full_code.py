def translate_dna_to_protein(dna_sequence):
    """
    Translates a DNA sequence into a protein sequence.

    The DNA sequence is read in codons (groups of three nucleotides), and each
    codon is mapped to its corresponding amino acid. The translation stops if
    a stop codon is encountered.

    The single-letter amino acid codes used in the translation map are:
        A: Alanine       G: Glycine       M: Methionine      S: Serine
        C: Cysteine      H: Histidine     N: Asparagine      T: Threonine
        D: Aspartic Acid I: Isoleucine    P: Proline         V: Valine
        E: Glutamic Acid K: Lysine        Q: Glutamine       W: Tryptophan
        F: Phenylalanine L: Leucine       R: Arginine        Y: Tyrosine

    Args:
        dna_sequence: A string representing the DNA sequence.

    Returns:
        A string representing the protein sequence.
    """

    # A dictionary that maps DNA codons to their corresponding amino acids.
    # '_' indicates a stop codon, which terminates translation.
    codon_to_amino_acid = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }

    protein_sequence = ""
    # Process the DNA sequence in chunks of three (codons).
    for i in range(0, len(dna_sequence), 3):
        # Extract a three-nucleotide codon from the DNA sequence.
        codon = dna_sequence[i:i+3]
        # Ensure the codon is of length 3 to avoid index errors.
        if len(codon) == 3:
            # Look up the amino acid for the codon and append it to the protein sequence.
            # If the codon is not found, an empty string is returned.
            amino_acid = codon_to_amino_acid.get(codon, '')
            if amino_acid == '_':
                # Stop translation if a stop codon is encountered.
                break
            protein_sequence += amino_acid
    return protein_sequence



def calculate_hamming_distance(slack_username, twitter_handle):
    """
    Calculates the Hamming distance between two strings.

    If the strings are of different lengths, the shorter string is padded
    with '#' to match the length of the longer string.

    Args:
        slack_username: The first string.
        twitter_handle: The second string.

    Returns:
        The Hamming distance between the two strings.
    """

    # Identify the lengths of both strings.
    length_slack = len(slack_username)
    length_twitter = len(twitter_handle)

    # Pad the shorter string to match the length of the longer string.
    if length_slack > length_twitter:
        twitter_handle = twitter_handle.ljust(length_slack, '#')
    elif length_twitter > length_slack:
        slack_username = slack_username.ljust(length_twitter, '#')

    hamming_distance = 0
    # Iterate through the strings and compare characters at each position.
    for i in range(len(slack_username)):
        if slack_username[i] != twitter_handle[i]:
            hamming_distance += 1

    return hamming_distance