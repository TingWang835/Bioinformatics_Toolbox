import sys
from Bio import SeqIO

# Experimental stiffness values (Persistence Length approximations)
# High values = Rigid (Promoter-like), Low values = Flexible (Histone-loving)
STIFFNESS_MAP = {
    "AAAAA": 5.0, "TTTTT": 5.0,  # Highly rigid poly-A/T tracts
    "GCGC": -3.0,   "CCGG": -3.0,    # Flexible, easily wrapped
    "DEFAULT": 0.0
}

def calculate_rigidity(sequence):
    scores = []
    window_size = 5
    for i in range(len(sequence) - window_size + 1):
        sub = str(sequence[i:i+window_size]).upper()
        # Find the best match in STIFFNESS_MAP or use default
        score = next((v for k, v in STIFFNESS_MAP.items() if k in sub), STIFFNESS_MAP["DEFAULT"])
        scores.append(score)
    return scores

def main(input_fasta, output_tsv):
    with open(output_tsv, "w") as f:
        f.write("id\tposition\trigidity_score\n")
        for record in SeqIO.parse(input_fasta, "fasta"):
            scores = calculate_rigidity(record.seq)
            for pos, score in enumerate(scores):
                f.write(f"{record.id}\t{pos}\t{score}\n")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])