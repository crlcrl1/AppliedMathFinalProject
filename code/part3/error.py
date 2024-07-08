from typing import List, Dict

from scipy import integrate as spi
from math import exp, sin, pi, cos
import matplotlib.pyplot as plt


def get_output(n: int, method: str):
    """
    Get the output of the solution from file
    """
    with open(f"output_{n}_{method}.txt", "r") as f:
        lines = f.readlines()
        temp = [list(map(float, line.strip().split())) for line in lines]
        for i in range(n - 1):
            temp[i].insert(0, 0.0)
            temp[i].append(0.0)
        temp.append([0.0] * (n + 1))
        temp.insert(0, [0.0] * (n + 1))
        return temp


def u_exact(x: float, y: float):
    """
    The exact solution
    """
    return exp(-2 * pi ** 2) * sin(pi * x) * sin(pi * y)


def u_exact_grad(x: float, y: float):
    """
    The gradient of the exact solution
    """
    return (pi * exp(-2 * pi ** 2) * sin(pi * y) * cos(pi * x),
            pi * exp(-2 * pi ** 2) * sin(pi * x) * cos(pi * y))


def calculate_error(n: int, method: str):
    """
    Calculate the error of the solution
    """
    u = get_output(n, method)
    error1 = 0.0
    error2 = 0.0
    h = 1 / n
    for i in range(n):
        for j in range(n):
            def f1(x: float, y: float):
                lt = u[i][j] * ((i + 1) * h - x) * ((j + 1) * h - y) / h ** 2
                lb = u[i + 1][j] * (x - i * h) * ((j + 1) * h - y) / h ** 2
                rt = u[i][j + 1] * ((i + 1) * h - x) * (y - j * h) / h ** 2
                rb = u[i + 1][j + 1] * (x - i * h) * (y - j * h) / h ** 2
                return (u_exact(x, y) - lt - lb - rt - rb) ** 2
            error1 += spi.dblquad(f1, j * h, (j + 1) * h, i * h, (i + 1) * h)[0]
            
            def f2(x: float, y: float):
                x_grad, y_grad = u_exact_grad(x, y)
                lt_x, lt_y = -u[i][j] * ((j + 1) * h - y) / h ** 2, -u[i][j] * ((i + 1) * h - x) / h ** 2
                lb_x, lb_y = u[i + 1][j] * ((j + 1) * h - y) / h ** 2, -u[i + 1][j] * (x - i * h) / h ** 2
                rt_x, rt_y = -u[i][j + 1] * (y - j * h) / h ** 2, u[i][j + 1] * ((i + 1) * h - x) / h ** 2
                rb_x, rb_y = u[i + 1][j + 1] * (y - j * h) / h ** 2, u[i + 1][j + 1] * (x - i * h) / h ** 2
                return (x_grad - lt_x - lb_x - rt_x - rb_x) ** 2 + (y_grad - lt_y - lb_y - rt_y - rb_y) ** 2
            error2 += spi.dblquad(f2, j * h, (j + 1) * h, i * h, (i + 1) * h)[0]
    return error1 ** 0.5, error2 ** 0.5


def write_to_file(n_list: List, res: Dict, filename: str):
    """
    Write the result to file
    """
    with open(filename, "w") as f:
        for method in res.keys():
            f.write(f"{method}:\n")
            for n, error in zip(n_list, res[method]):
                f.write(f"n = {n}\terror = {error}\n")
            f.write("\n")


def plot(n_list: List, errors: Dict):
    """
    Plot the error
    """
    for method in errors.keys():
        plt.loglog(n_list, errors[method], marker="o", label=method)
    plt.xlabel("n")
    plt.ylabel("error")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    n_list = [4, 8, 16, 32, 64, 128]
    method_list = ["explicit", "Crank-Nicolson", "implicit"]
    errors1_dict = {}
    errors2_dict = {}
    for method in method_list:
        errors1 = []
        errors2 = []
        for n in n_list:
            error1, error2 = calculate_error(n, method)
            errors1.append(error1)
            errors2.append(error2)
        errors1_dict[method] = errors1
        errors2_dict[method] = errors2
    
    write_to_file(n_list, errors1_dict, "error_L2.txt")
    write_to_file(n_list, errors2_dict, "error_H1.txt")
    
    plot(n_list, errors1_dict)
    plot(n_list, errors2_dict)
