import numpy as np
from hmmlearn import hmm

# Parameters
n_states = 2  # Fair, Unfair
n_observations = 6  # Dice faces 1-6
length = 20  # Sequence length

# State names
states = ["Fair", "Unfair"]

# Transition matrix: probability of switching between fair/unfair
transition_matrix = np.array([
    [0.95, 0.05],  # Fair -> Fair, Fair -> Unfair
    [0.10, 0.90]   # Unfair -> Fair, Unfair -> Unfair
])

# Emission matrix: probability of each dice face in each state
emission_matrix = np.array([
    [1/6]*6,                # Fair dice: uniform probability
    [0.05, 0.05, 0.05, 0.05, 0.05, 0.75]  # Unfair dice: 6 is much more likely
])

# Initial state probabilities
start_prob = np.array([0.5, 0.5])

# Create HMM model
model = hmm.MultinomialHMM(n_components=n_states, init_params="")
model.startprob_ = start_prob
model.transmat_ = transition_matrix
model.emissionprob_ = emission_matrix

# Generate a random sequence
X, Z = model.sample(length)
sequence = X.flatten() + 1  # Dice faces (1-6)

print("Generated sequence:", sequence)
print("True states:", [states[z] for z in Z])

# Use Viterbi to decode the most likely state sequence
logprob, decoded_states = model.decode(X, algorithm="viterbi")
print("Viterbi decoded states:", [states[s] for s in decoded_states])

# Show matrices
print("\nTransition matrix:\n", transition_matrix)
print("\nEmission matrix:\n", emission_matrix)