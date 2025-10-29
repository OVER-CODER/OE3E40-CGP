# phylo_tree_upgma_modern_fixed_v3.py

from Bio import Entrez, SeqIO, Phylo
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
import ssl, certifi

# ---- SSL & NCBI setup ----
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())
Entrez.email = "your_email@example.com"  # replace with your email

# ---- GenBank data ----
data = [
    ("German_Neanderthal", "AF011222"),
    ("Russian_Neanderthal", "AF254446"),
    ("European_Human", "X90314"),
    ("Mountain_Gorilla_Rwanda", "AF089820"),
    ("Chimp_Troglodytes", "AF176766"),
]

# ---- Fetch sequences ----
records = []
for header, acc in data:
    print(f"Fetching {header} ({acc}) ...")
    handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
    rec = SeqIO.read(handle, "fasta")
    handle.close()
    rec.id = header
    rec.name = header
    rec.description = header
    records.append(rec)

# ---- Pairwise aligner ----
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.open_gap_score = -2
aligner.extend_gap_score = -0.5

def pairwise_align_to_gapped(seq1, seq2):
    """Perform a global alignment and return two equal-length sequences."""
    aln = aligner.align(seq1, seq2)[0]
    aligned_seq1 = aln.target
    aligned_seq2 = aln.query
    return str(aligned_seq1), str(aligned_seq2)

# ---- Progressive alignment ----
aligned_records = [records[0]]

for rec in records[1:]:
    ref_seq = str(aligned_records[0].seq)
    new_seq = str(rec.seq)
    ref_aln, new_aln = pairwise_align_to_gapped(ref_seq, new_seq)

    # Pad previous sequences to match new alignment length
    diff = len(ref_aln) - len(ref_seq)
    if diff > 0:
        for i, r in enumerate(aligned_records):
            aligned_records[i].seq = Seq(str(r.seq) + "-" * diff)

    aligned_records[0].seq = Seq(ref_aln)
    rec.seq = Seq(new_aln)
    aligned_records.append(rec)

# ---- Normalize lengths ----
max_len = max(len(r.seq) for r in aligned_records)
for r in aligned_records:
    if len(r.seq) < max_len:
        r.seq = Seq(str(r.seq) + "-" * (max_len - len(r.seq)))

alignment = MultipleSeqAlignment(aligned_records)
print(f"\n✅ Alignment complete — all sequences length = {max_len}")

# ---- Distance computation ----
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

# ---- Jukes–Cantor correction (safe version) ----
def jukes_cantor(p):
    p = min(max(p, 0.0), 0.7499)  # clamp to avoid log domain errors
    return -3/4 * np.log(1 - 4*p/3)

for i in range(len(dm.names)):
    for j in range(len(dm.names)):
        if i != j:
            dm[i, j] = jukes_cantor(dm[i, j])

print("\nJukes–Cantor distance matrix:")
print(dm)

# Replace NaN or Inf with max finite value
finite_vals = [v for row in dm for v in row if np.isfinite(v)]
max_finite = max(finite_vals) if finite_vals else 1.0

for i in range(len(dm.names)):
    for j in range(len(dm.names)):
        if not np.isfinite(dm[i, j]):
            dm[i, j] = max_finite

# ---- Tree construction (UPGMA) ----
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)

# ---- Visualization ----
print("\n🧬 ASCII Phylogenetic Tree:")
Phylo.draw_ascii(tree)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=ax, do_show=False)
plt.title("UPGMA Phylogenetic Tree (Jukes–Cantor)")
plt.ylabel("Evolutionary Distance")
plt.show()
