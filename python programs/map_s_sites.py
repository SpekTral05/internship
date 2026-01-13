import argparse

def parse_fasta(fasta_file):
    """Returns dict: accession -> full protein sequence"""
    sequences = {}
    current_acc = None
    seq_chunks = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_acc:
                    sequences[current_acc] = "".join(seq_chunks)
                current_acc = line.split("|")[1] if "|" in line else line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if current_acc:
            sequences[current_acc] = "".join(seq_chunks)

    return sequences


def extract_s_number(mod_string):
    """
    Extracts the number after 'S' (e.g. S10 → 10)
    Returns None if not found
    """
    s_index = mod_string.find("S")
    if s_index == -1:
        return None

    num = ""
    i = s_index + 1
    while i < len(mod_string) and mod_string[i].isdigit():
        num += mod_string[i]
        i += 1

    if num == "":
        return None

    return int(num)


def find_absolute_s_position(peptide, s_number, protein_sequence):
    pep_start = protein_sequence.find(peptide)
    if pep_start == -1:
        return None

    return pep_start + s_number


def process_input(input_file, fasta_file):
    proteins = parse_fasta(fasta_file)

    with open(input_file, "r") as f:
        next(f)  # skip header
        for line in f:
            cols = line.strip().split("\t")
            peptide = cols[0]
            mods = cols[1]
            accession = cols[2]

            s_number = extract_s_number(mods)
            if s_number is None:
                continue

            protein_seq = proteins.get(accession)
            if not protein_seq:
                print(f"{accession}: not found in FASTA")
                continue

            abs_pos = find_absolute_s_position(peptide, s_number, protein_seq)

            if abs_pos is None:
                print(f"{accession}: peptide not found")
            else:
                print(f"{accession}\tS{abs_pos}")


def main():
    parser = argparse.ArgumentParser(description="Map peptide S-sites to protein positions (no regex)")
    parser.add_argument("input", help="Input TSV file")
    parser.add_argument("fasta", help="FASTA file")
    args = parser.parse_args()

    process_input(args.input, args.fasta)


if __name__ == "__main__":
    main()
