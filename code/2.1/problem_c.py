from shared import *

from numpy.typing import NDArray
import numpy as np
import matplotlib.pyplot as plt


def plot(res: NDArray):
    n = len(res) + 1
    x = np.linspace(0, 1, n + 1)
    plt.plot(x, np.concatenate(([0], res, [0])), label=f"n={n}")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    n_list = [10, 20, 40, 80]
    for n in n_list:
        res = solve(n)
        plot(res)
