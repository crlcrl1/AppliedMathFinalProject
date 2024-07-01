from shared import *


if __name__ == "__main__":
    n_list = [10, 20, 40, 80]
    for n in n_list:
        res = solve(n)
        plot(np.linspace(0, 1, n + 1), res)
