from Bio import SeqIO
from collections import Counter

def read_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    return records[0]

def calculate_gc_content(sequence):
    sequence = sequence.upper()
    g = sequence.count("G")
    c = sequence.count("C")
    return round(((g + c) / len(sequence)) * 100, 2)

def count_nucleotides(sequence):
    sequence = sequence.upper()
    return Counter(sequence)

def find_orfs(sequence, min_length=90):
    sequence = str(sequence).upper()
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    orfs = []

    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]

            if codon == start_codon:
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop = sequence[j:j+3]

                    if stop in stop_codons:
                        length = j + 3 - i
                        if length >= min_length:
                            orfs.append({
                                "start": i + 1,
                                "end": j + 3,
                                "length": length,
                                "frame": frame + 1,
                                "sequence": sequence[i:j+3]
                            })
                        break
    return orfs

def save_report(record, nucleotide_counts, gc_content, orfs):
    with open("results/genome_report.txt", "w") as report:
        report.write("Micobial Genome Explorer Report\n")
        report.write("================================\n\n")
        report.write(f"Sequence ID: {record.id}\n")
        report.write(f"Genome length: {len(record.seq)} base pairs\n")
        report.write(f"GC content: {gc_content}%\n\n")

        report.write("Nucleotide Counts:\n")
        for base, count in nucleotide_counts.items():
            report.write(f"{base}: {count}\n")

        report.write(f"\nOpen Reading Frames found: {len(orfs)}\n\n")
        for index, orf in enumerate(orfs[:20], start=1):
            report.write(f"ORF {index}: Start {orf['start']}, End {orf['end']}, Length {orf['length']}, Frame {orf['frame']}\n")

def main():
    fasta_file = "data/ecoli_k12.fasta"

    record = read_fasta(fasta_file)
    sequence = record.seq

    nucleotide_counts = count_nucleotides(sequence)
    gc_content = calculate_gc_content(sequence)
    orfs = find_orfs(sequence)

    print("Analysis complete!")
    print(f"Sequence ID: {record.id}")
    print(f"Genome length: {len(sequence)} base pairs")
    print(f"GC content: {gc_content}%")
    print(f"ORFs found: {len(orfs)}")

    save_report(record, nucleotide_counts, gc_content, orfs)
    print("Report saved to results/genome_report.txt")

if __name__ == "__main__":
    main()
