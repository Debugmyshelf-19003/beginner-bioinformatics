def read_fasta(file_path):
    sequence = ""
    with open(file_path, "r") as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip().upper()
    return sequence


def gc_content(sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return round((gc_count / len(sequence)) * 100, 2)


def nucleotide_counts(sequence):
    return {
        "A": sequence.count("A"),
        "T": sequence.count("T"),
        "G": sequence.count("G"),
        "C": sequence.count("C")
    }


def reverse_complement(sequence):
    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    return "".join(complement[base] for base in reversed(sequence))


sequence = read_fasta("sample_sequence.fasta")

print("Sequence length:", len(sequence))
print("GC content:", gc_content(sequence), "%")
print("Nucleotide counts:", nucleotide_counts(sequence))
print("Reverse complement:", reverse_complement(sequence))
