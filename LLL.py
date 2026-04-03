# LLL Algorithm Implementation following An_Introduction_to_Mathematical_Cryptography by Silverman, Figure 7.8

import numpy as np

# The Gram-Schmidt process is used to orthogonalize the basis vectors, which is essential for the LLL algorithm to check the Lovász condition and perform size reduction effectively.
# It must not return the orthonormal basis, but rather the orthogonalized basis, as the LLL algorithm relies on the lengths of these vectors to determine when to swap and reduce.
def gram_schmidt(B):
    B_star = np.zeros_like(B, dtype=float)
    B_star[0] = B[0]
    for i in range(1, len(B)):
        B_star[i] = B[i]
        for j in range(i):
            mu = np.dot(B[i], B_star[j]) / np.dot(B_star[j], B_star[j])
            B_star[i] -= mu * B_star[j]
    return B_star


# The LLL algorithm takes a basis B of a lattice and a parameter delta (which controls the "reduction" level, typically set to 0.75) and returns a reduced basis of the same lattice.
def LLL(B, delta=0.75):
    B_star = gram_schmidt(B)
    k = 1
    n = len(B)
    while k < n:
        for j in range(k-1, -1, -1):
            mu_kj = np.dot(B[k], B_star[j]) / np.dot(B_star[j], B_star[j]) # Projection of B[k] onto B_star[j], normalized by the length of B_star[j], as indicated in the foot note.
            B[k] -= round(mu_kj) * B[j] # Size Reduction
                
        # At this poing, we have reduced B[k] with respect to all previous vectors. Now we check the Lovász condition to decide if we need to swap B[k] and B[k-1].
        if np.dot(B_star[k], B_star[k]) >= (delta - ((np.dot(B[k], B_star[k-1]))/np.dot(B_star[k-1], B_star[k-1]))**2 ) * np.dot(B_star[k-1], B_star[k-1]):
            k += 1
        else:
            B[[k, k-1]] = B[[k-1, k]] # Swap B[k] and B[k-1]
            B_star = gram_schmidt(B)  
            k = max(k-1, 1)
    return B
            
if __name__ == "__main__":
    # Example usage
    M = np.array([
    [19,  2, 32, 46,  3, 33],
    [15, 42, 11,  0,  3, 24],
    [43, 15,  0, 24,  4, 16],
    [20, 44, 44,  0, 18, 15],
    [ 0, 48, 35, 16, 31, 31],
    [48, 33, 32,  9,  1, 29]]) # Silverman Example 7.75
    
    print("LLL Reduction:")
    LM = LLL(M)
    print(LM)
    
    print(np.linalg.det(M) == np.linalg.det(LM)) # The determinant (volume of the fundamental parallelepiped) should be the same as it's an invariant of the lattice.