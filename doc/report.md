# 应用数学导论大作业报告

## 一维有限元方法

### 有限元方法

(a) 我们有：
$$
\int_0^L \frac{{\rm d}u}{{\rm d}x} \frac{{\rm d}v}{{\rm d}x} {{\rm d}x} =
-\int_0^L \frac{{\rm d}^2u}{{\rm d}x^2}v {{\rm d}x}
$$
从而
$$
\int_0^L \left(\frac{{\rm d}^2u}{{\rm d}x^2} + f\right)v {{\rm d}x}
$$
取
$$
v(x)=
\begin{cases}
1,& t\le x\le t+\epsilon\\
0,& else 
\end{cases}
$$
其中 $0<t<1$，有
$$
\int_0^L \left(\frac{{\rm d}^2u}{{\rm d}x^2} + f\right)v {{\rm d}x} = \epsilon \left(\frac{{\rm d}^2u}{{\rm d}x^2} + f\right)(\xi) = 0
$$
其中 $t\le\xi\le t + \epsilon$，令 $\epsilon \to 0$，有
$$
\left(\frac{{\rm d}^2u}{{\rm d}x^2} + f\right)(t) = 0
$$
由 $t$ 的任意性，
$$
\frac{{\rm d}^2u}{{\rm d}x^2} + f\equiv0,\,\,\forall x\in(0,L)
$$
（b）由分部积分：
$$
\int_{\Omega}\frac{{\rm d}(u-u_h)}{{\rm d}x} \frac{{\rm d}v_h}{{\rm d}x}{\rm d}x = 
(u-u_h)\frac{{\rm d}v_h}{{\rm d}x}|_0^L - \int_{\Omega}(u-u_h)\frac{{\rm d}^2v_h}{{\rm d}x^2}{\rm d}x
$$
由于 $(u-u_h)(0) = (u-u_h)(L)=0$ 和 $\frac{{\rm d}^2v_h}{{\rm d}x^2}=0$，有原式为 0

（c）问题（1.1）的变分形式为：
$$
\int_0^1\frac{{\rm d}u}{{\rm d}x} \frac{{\rm d}v}{{\rm d}x}{\rm d}x = \int_0^1u{\rm d}x
$$
单元个数为 5 时
$$
A = 
\begin{pmatrix}
\frac{1}{h_1}+\frac{1}{h_2} & -\frac{1}{h_2} \\
-\frac{1}{h_2} & \frac{1}{h_2}+\frac{1}{h_3} & -\frac{1}{h_3}\\
& -\frac{1}{h_3} & \frac{1}{h_3}+\frac{1}{h_4} & -\frac{1}{h_4}\\
& & -\frac{1}{h_4} & \frac{1}{h_4}+\frac{1}{h_5} \\
\end{pmatrix}
$$
