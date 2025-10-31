import numpy as np
import os
import time
from Bio import SeqIO
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

# ==============================================================================
# 1. EUCLIDEAN DISTANCE CALCULATION
# ==============================================================================

def getEDistance(A, B):
    """
    Get Euclidean distance between corresponding rows of matrices A and B.

    Rows represent data points; columns represent dimensions. A and B must
    have the same number of columns (dimensions).

    Args:
        A (np.ndarray): The first matrix of data points.
        B (np.ndarray): The second matrix of data points.

    Returns:
        np.ndarray: A 1D array where each element is the Euclidean distance
                    between the corresponding row in A and B.
    """
    # Euclidean distance = sqrt(sum((A - B)**2, axis=1))
    return np.sqrt(np.sum((A - B)**2, axis=1))

# ==============================================================================
# 2. MOMENT VECTOR CALCULATION (Power Spectrum Method)
# ==============================================================================

def GetMomentVectorPS(seq):
    """
    Calculates a moment vector based on the power spectrum (PS) 
    of nucleotide indicator sequences.

    Args:
        seq (str): The input DNA sequence (e.g., 'ATGCAGTC').

    Returns:
        np.ndarray: A 1D array (vector) containing the first three moments 
                    for A, C, G, and T, in the order [M1A, M1C, M1G, M1T, 
                    M2A, M2C, M2G, M2T, M3A, M3C, M3G, M3T].
    """
    n = len(seq)
    seq_upper = seq.upper()
    
    # Initialize indicator sequences and counts
    uA, uC, uG, uT = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    nA, nC, nG, nT = 0, 0, 0, 0
    
    # Get indicator sequences and counts
    for i in range(n):
        nu = seq_upper[i]
        if nu == 'A':
            uA[i] = 1; nA += 1
        elif nu == 'C':
            uC[i] = 1; nC += 1
        elif nu == 'G':
            uG[i] = 1; nG += 1
        elif nu == 'T':
            uT[i] = 1; nT += 1
        
    # Discrete Fourier Transforms (DFT/FFT)
    UA, UC, UG, UT = np.fft.fft(uA), np.fft.fft(uC), np.fft.fft(uG), np.fft.fft(uT)
    
    # Exclude the first term (DC component, index 0)
    UA, UC, UG, UT = UA[1:], UC[1:], UG[1:], UT[1:]
    
    # Power spectrums (PS = |DFT|^2)
    PSA = np.abs(UA) ** 2     
    PSC = np.abs(UC) ** 2     
    PSG = np.abs(UG) ** 2     
    PST = np.abs(UT) ** 2     
    
    # Initialize moment vectors
    MA, MC, MG, MT = np.zeros(3), np.zeros(3), np.zeros(3), np.zeros(3)
    
    PS_list = [PSA, PSC, PSG, PST]
    n_list = [nA, nC, nG, nT]
    M_list = [MA, MC, MG, MT]

    # Calculate moments (j=1 to 3)
    for ps_vec, count, moment_vec in zip(PS_list, n_list, M_list):
        if count > 0 and n - count > 0:
            for j in range(1, 4):  # j goes from 1 to 3 (moments M1, M2, M3)
                sum_ps_j = np.sum(ps_vec ** j)
                num_factor = count * (n - count)
                den_factor = (count ** j) * ((n - count) ** j)
                
                # Store moment at index j-1
                moment_vec[j - 1] = num_factor * sum_ps_j / den_factor
        
    # Moment vector assembly: M1A, M1C, M1G, M1T, M2A, ..., M3T
    v = np.concatenate([
        MA[[0]], MC[[0]], MG[[0]], MT[[0]],
        MA[[1]], MC[[1]], MG[[1]], MT[[1]],
        MA[[2]], MC[[2]], MG[[2]], MT[[2]]
    ])
    
    return v

# ==============================================================================
# 3. UPGMA PHYLOGENETIC ANALYSIS
# ==============================================================================

def TestUPGMA(name):
    """
    Reads a FASTA file, calculates normalized moment vectors, computes a distance 
    matrix, and builds and plots a UPGMA phylogenetic tree.

    Args:
        name (str): The base name of the FASTA file (e.g., 'sequences' 
                    for 'sequences.fasta').
    """
    
    start_time = time.time()
    file = f"{name}.fasta"
    
    # Read FASTA file
    try:
        sequences = list(SeqIO.parse(file, "fasta"))
    except FileNotFoundError:
        print(f"\nError: FASTA file '{file}' not found. Please ensure it is in the current directory.")
        return

    len_seqs = len(sequences)
    if len_seqs < 2:
        print(f"Error: Not enough sequences ({len_seqs}) in '{file}' to build a tree.")
        return

    # Analyze sequence lengths
    len_x = [len(s.seq) for s in sequences]
    min_len, max_len = min(len_x), max(len_x)
    print(f'\nMin sequence length: {min_len}')
    print(f'Max sequence length: {max_len}')
    
    # Get and normalize moment vectors
    b = max_len
    v = []
    print('Calculating Moment Vectors...')
    for seq_record in sequences:
        moment_vec = GetMomentVectorPS(str(seq_record.seq))
        # Normalize by max length: 1/b * MomentVector
        v.append(moment_vec / b) 
        
    V_matrix = np.array(v)
    
    # Calculate Condensed Euclidean Distance Matrix (d)
    print('Calculating Distance Matrix...')
    # pdist computes D(i, j) for i > j (lower triangle as a vector)
    d = pdist(V_matrix, metric='euclidean')
    
    # Phylogenetic tree construction (UPGMA Linkage)
    print('Building Phylogenetic Tree (UPGMA)...')
    Z = linkage(d, method='average')
    
    # Plotting the Dendrogram
    labels = [s.id for s in sequences]
    
    plt.figure(figsize=(10, len(labels) * 0.3 + 2)) 
    
    dendrogram(
        Z,
        orientation='left',  # Matches MATLAB's 'orient', 'left'
        labels=labels,
        leaf_font_size=8
    )
    
    plt.title('Similarity distance using our new method with UPGMA', 
              fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # Timer end
    t_end = time.time() - start_time
    minutes = int(t_end // 60)
    seconds = t_end % 60
    
    print(f'\nAnalysis completed in: {minutes} minutes and {seconds:.4f} seconds')


# ==============================================================================
# 4. MAIN INTERACTIVE MENU
# ==============================================================================

def main_menu():
    """Presents a text-based menu for selecting a dataset and runs the analysis."""
    
    # Clear console ('clc' equivalent)
    os.system('cls' if os.name == 'nt' else 'clear') 
    
    print("Welcome to the UPGMA Phylogenetic Tree Builder!")

    condition = True
    while condition:
        
        print("\n" + "="*35)
        print("Please choose a data set to analyze:")
        print("  1. Mammals (Mammals.fasta)")
        print("  2. Influenza A virus (Influenza.fasta)")
        print("  3. Coronavirus (Corona.fasta)")
        print("  4. HRV (HRV.fasta)")
        print("  5. Bacteria (Bacteria.fasta)")
        print("  6. Quit")
        print("="*35)
        
        try:
            choice = input("Enter your choice (1-6): ").strip()
            choice_int = int(choice)
            
            dataset_map = {
                1: 'Mammals', 2: 'Influenza', 3: 'Corona', 
                4: 'HRV', 5: 'Bacteria'
            }
            
            if choice_int in dataset_map:
                TestUPGMA(dataset_map[choice_int])
                
            elif choice_int == 6:
                print('\nYou closed the menu!')
                condition = False
            else:
                print("Invalid choice. Please enter a number between 1 and 6.")
                
        except ValueError:
            print("Invalid input. Please enter a valid number.")
        
    print('\nExiting program.')

if __name__ == "__main__":
    main_menu()