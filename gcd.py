import LLL
import numpy as np
import math

def gcd(numbers):
    
    basis = np.zeros((len(numbers), len(numbers)+1))
    for i in range(0, len(numbers)):
        row = np.zeros(len(numbers)+1)
        row[i] = 1
        row[len(numbers)] = numbers[i]
        basis[i] = row
    
    reduced_basis = LLL.LLL(basis)
    #print(reduced_basis)
    
    for i in range(len(numbers)):
        if reduced_basis[i][len(numbers)] == 0:
            continue
        
        sum = 0
        for j in range(len(numbers)+1):
            if j != len(numbers):
               sum += reduced_basis[i][j] * numbers[j]
            else:
                gcd = reduced_basis[i][j]
                if sum ==  gcd:
                    return reduced_basis[i] if gcd > 0 else -1 * reduced_basis[i]
    return [0]*(len(numbers)+1)

if __name__ == "__main__":
    arr = [48, 18, 30, 12, 120]
    lincomb = gcd(arr)
    
    for i in range(len(lincomb)):
        if i != len(arr):
            print(f"{arr[i]} * ({lincomb[i]}) + ", end="")
        else:
            print(f"= {lincomb[i]}")