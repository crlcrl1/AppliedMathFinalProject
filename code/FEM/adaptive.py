import scipy.integrate as spi

from shared import *


node_list: List[float] = np.linspace(0, 1, 11).tolist()
cache: List[float] = []


def f(x):
    return np.exp(-100 * (x - 0.5) ** 2)


def setup():
    global node_list
    global cache
    node_list = np.linspace(0, 1, 11).tolist()
    cache = []


def calculate_integral():
    n = len(node_list) - 1
    cache.clear()
    for i in range(n):
        cache.append(spi.quad(lambda x: f(x) ** 2, node_list[i], node_list[i + 1])[0])


def densify(alpha: float):
    global node_list
    max_v = max(cache)
    new_node_list = []
    for i, v in enumerate(cache):
        new_node_list.append(node_list[i])
        if v > alpha * max_v:
            new_node_list.append((node_list[i] + node_list[i + 1]) / 2)
    new_node_list.append(node_list[-1])
    node_list = new_node_list


def generate_A(node_list: List[float]):
    n = len(node_list) - 1
    temp = np.zeros((n - 1, n - 1), dtype=np.float64)
    h_list = [1 / (node_list[i + 1] - node_list[i]) for i in range(n)]
    for i in range(n - 1):
        temp[i, i] = h_list[i] + h_list[i + 1]
        if i < n - 2:
            temp[i, i + 1] = -h_list[i + 1]
            temp[i + 1, i] = -h_list[i + 1]
    return temp


def generate_b(node_list: List[float]):
    n = len(node_list) - 1
    h_list = [node_list[i + 1] - node_list[i] for i in range(n)]
    return np.array([f(node_list[i + 1]) * (h_list[i] + h_list[i + 1]) / 2 for i in range(n - 1)],
                    dtype=np.float64)


def solve(n: int, alpha: float):
    setup()
    while len(node_list) <= n:
        calculate_integral()
        densify(alpha)
    A = generate_A(node_list)
    b = generate_b(node_list)
    plot(node_list, chase_method(A, b))


if __name__ == "__main__":
    solve(80, 0.1)
