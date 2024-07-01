from typing import List
import matplotlib.pyplot as plt

from shared import *


def L2_error(u: NDArray) -> float:
    n = len(u) + 1
    u = np.concatenate(([0], u, [0]))
    error = 0.0
    for i in range(n):
        a1 = i / n
        a2 = (i + 1) / n
        a = 1 / 2 - (u[i + 1] - u[i]) * n
        error += a ** 2 * (a2 - a1) - a * (a2 ** 2 - a1 ** 2) + 1 / 3 * (a2 ** 3 - a1 ** 3)
    return np.sqrt(error)


def make_table(n_list: List, errors: List):
    """
    Print a table in the format of markdown
    """
    print("|  n  | L2 error |")
    print("|-----|----------|")
    for n, error in zip(n_list, errors):
        print(f"| {n:^3} | {error:.6f} |")


if __name__ == "__main__":
    n_list = [2, 4, 8, 16, 32, 64, 128]
    errors = []
    for n in n_list:
        u = solve(n)
        errors.append(L2_error(u))
    make_table(n_list, errors)
    # plot
    plt.loglog(n_list, errors, marker="o")
    plt.xlabel("n")
    plt.ylabel("L2 error")
    plt.show()
