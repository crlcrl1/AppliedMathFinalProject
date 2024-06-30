"""
This module contains shared functions for the 2.1 task.
"""

import numpy as np

def generate_A(n: int):
    temp = np.zeros((n - 1, n - 1), dtype=np.float64)
    for i in range(n - 1):
        temp[i, i] = 2 * n
        if i < n - 2:
            temp[i, i + 1] = -n
            temp[i + 1, i] = -n
    return temp


def generate_b(n: int):
    return np.array([1 / n] * (n - 1), dtype=np.float64)


def solve(n: int):
    A = generate_A(n)
    b = generate_b(n)
    for i in range(1, n - 1):
        factor = A[i, i - 1] / A[i - 1, i - 1]
        A[i, i] -= factor * A[i - 1, i]
        b[i] -= factor * b[i - 1]
    for i in range(n - 3, -1, -1):
        factor = A[i, i + 1] / A[i + 1, i + 1]
        b[i] -= factor * b[i + 1]
    return b / np.diag(A)
