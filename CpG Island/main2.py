import numpy as np
from hmmlearn import hmm
import matplotlib.pyplot as plt
import os

# Read the file
def read_file(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
    return ''.join(content.split())

def cpgisland(sequence, plot=False):
    """
    Simplified CpG island detection
    Returns start and stop positions of potential CpG islands
    """
    window_size = 200
    min_gc_content = 0.5
    min_oe_ratio = 0.6
    
    sequence = sequence.upper()
    n = len(sequence)
    gc_content = []
    oe_ratio = []
    
    for i in range(n - window_size + 1):
        window = sequence[i:i + window_size]
        
        # Calculate GC content
        gc = (window.count('G') + window.count('C')) / window_size
        gc_content.append(gc)
        
        # Calculate observed/expected CpG ratio
        cpg_observed = window.count('CG')
        c_count = window.count('C')
        g_count = window.count('G')
        cpg_expected = (c_count * g_count) / window_size if window_size > 0 else 0
        ratio = cpg_observed / cpg_expected if cpg_expected > 0 else 0
        oe_ratio.append(ratio)
    
    # Find regions that meet criteria
    islands = []
    start = None
    
    for i in range(len(gc_content)):
        if gc_content[i] >= min_gc_content and oe_ratio[i] >= min_oe_ratio:
            if start is None:
                start = i
        elif start is not None:
            islands.append((start, i))
            start = None
    
    if start is not None:
        islands.append((start, len(gc_content) - 1))
    
    if plot:
        plt.figure(figsize=(15, 6))
        plt.subplot(2, 1, 1)
        plt.plot(gc_content)
        plt.axhline(y=min_gc_content, color='r', linestyle='--')
        plt.title('GC Content')
        plt.ylabel('GC Content')
        
        plt.subplot(2, 1, 2)
        plt.plot(oe_ratio)
        plt.axhline(y=min_oe_ratio, color='r', linestyle='--')
        plt.title('Observed/Expected CpG Ratio')
        plt.xlabel('Position')
        plt.ylabel('O/E Ratio')
        
        plt.tight_layout()
        plt.show()
    
    return islands

def main():
    # Read sequence data
    file_path = 'hmr195.fa'
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        print("Please make sure hmr195.fa exists in the current directory")
        return
    g = read_file(file_path)
    
    # Get CpG islands
    start_stop = cpgisland(g)
    cpgisland(g, plot=True)
    
    # Convert sequence to numeric codes
    seq_mapping = {'A': 1, 'C': 2, 'G': 3, 'T': 4}
    SEQ = np.array([seq_mapping.get(base, 0) for base in g[:-1]])
    
    # Create STAT array with different regions
    region_mappings = [
        (0, 695, {'A': 5, 'C': 6, 'G': 7, 'T': 8}),
        (695, 1031, {'A': 1, 'C': 2, 'G': 3, 'T': 4}),
        (1031, 1054, {'A': 5, 'C': 6, 'G': 7, 'T': 8}),
        (1054, 1291, {'A': 1, 'C': 2, 'G': 3, 'T': 4}),
        (1291, 1383, {'A': 5, 'C': 6, 'G': 7, 'T': 8}),
        (1383, 1612, {'A': 1, 'C': 2, 'G': 3, 'T': 4}),
        (1612, len(g) - 1, {'A': 5, 'C': 6, 'G': 7, 'T': 8})
    ]
    
    STAT = np.zeros(len(g) - 1)
    for start, end, mapping in region_mappings:
        for i in range(start, end):
            if i < len(g) - 1:
                STAT[i] = mapping.get(g[i], 0)
    
    # Hidden Markov Model analysis
    # Convert sequences to counts format
    n_features = 4  # A, C, G, T
    n_samples = len(SEQ)
    
    # Convert SEQ to counts format (window-based approach)
    window_size = 1  # Using single nucleotides
    n_windows = n_samples - window_size + 1
    X = np.zeros((n_windows, n_features), dtype=int)
    
    for i in range(n_windows):
        window = SEQ[i:i + window_size]
        for nucleotide in window:
            if 1 <= nucleotide <= 4:
                X[i, int(nucleotide-1)] += 1
    
    # Initialize and train HMM
    model = hmm.MultinomialHMM(n_components=2, n_iter=100)
    
    # Fit the model
    model.fit(X)
    
    # Get transition and emission matrices
    TRANS = model.transmat_
    EMIS = model.emissionprob_
    
    # Find most likely state sequence
    likelystates = model.predict(X)
    
    # Get state probabilities
    PSTATES = model.predict_proba(X)
    
    print("Transition Matrix:")
    print(TRANS)
    print("\nEmission Matrix:")
    print(EMIS)
    print("\nFirst few likely states:")
    print(likelystates[:10])
    print("\nFirst few state probabilities:")
    print(PSTATES[:5])

if __name__ == "__main__":
    main()
