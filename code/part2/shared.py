"""
This module contains shared functions for the 2.1 task.
"""
from typing import List

import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt


def chase_method(A: NDArray, b: NDArray) -> NDArray:
    """
    Use chase method to solve a liner equation
    """
    n = len(b) + 1
    for i in range(1, n - 1):
        factor = A[i, i - 1] / A[i - 1, i - 1]
        A[i, i] -= factor * A[i - 1, i]
        b[i] -= factor * b[i - 1]
    for i in range(n - 3, -1, -1):
        factor = A[i, i + 1] / A[i + 1, i + 1]
        b[i] -= factor * b[i + 1]
    return b / np.diag(A)


def generate_A(n: int):
    """
    Generate matrix A for the equidistant partition
    """
    temp = np.zeros((n - 1, n - 1), dtype=np.float64)
    for i in range(n - 1):
        temp[i, i] = 2 * n
        if i < n - 2:
            temp[i, i + 1] = -n
            temp[i + 1, i] = -n
    return temp


def generate_b(n: int):
    """
    Generate vector b for the equidistant partition
    """
    return np.array([1 / n] * (n - 1), dtype=np.float64)


def solve(n: int):
    """
    Solve the equation with n intervals and equidistant partition
    """
    A = generate_A(n)
    b = generate_b(n)
    return chase_method(A, b)


def plot(x: NDArray | List, y: NDArray):
    """
    Plot the result
    """
    n = len(y) + 1
    y = np.concatenate(([0], y, [0]))
    plt.plot(x, y, label=f"n={n}")
    plt.scatter(x, y, s=5)
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.show()
