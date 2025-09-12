import matplotlib.pyplot as plt


seq = "GTGCGTATGACTGCAGAAATTTTGGTAGTTGGTTCAATTCTGCCAAGTGACAAATTTTTTTGTTTTTTGTGAGTA"

mapping = {
    "A": 0.1260,
    "T": 0.1335,
    "G": 0.0806,
    "C": 0.1340
}

res = []

for s in seq:
    if s in mapping:
        res.append(mapping[s])

plt.plot(res)
plt.show()

print(res)