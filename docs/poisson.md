# What the moving window measures

Run the following from the repository root.

~~~sh
python -m primes_to_zeros poisson 64 128 256 --relative-width 0.125 --modes 64
~~~

Let $c$ be the center on the integer axis, let $r$ be the relative width, and set $h=rc$. The window uses the centered cubic B-spline $B$, the convolution of four copies of the unit-mass interval $[-1/2,1/2]$. It is even, has integral one, has two continuous derivatives, and vanishes outside $[-2,2]$. In log coordinates the test function is

```math
\phi_c(u)=\sqrt{e^u/c}\,B\big((e^u-c)/h\big).
```

With $r$ fixed, this is a fixed shape translated by $\log c$. The experiment requires $c-2h\ge1$ to keep the full window on the positive half of the log axis. The weight of each spike times its window value simplifies to $\Lambda(n)B((n-c)/h)/\sqrt c$. [Unser, Aldroubi and Eden, *B-Spline Signal Processing*, Part I](https://bigwww.epfl.ch/publications/unser9301.pdf)

At integer centre $c=p^j$ and absolute width $h=1/4$, the spline support is $[c-1/2,c+1/2]$. Its only integer is $p^j$. Since $B(0)=2/3$, the prime-side sum is $2\log p/(3\sqrt{p^j})$.

Weight the same spline by the primitive character coefficient $u_j(p)=\sum_\chi\chi(p)^j$ from the explicit formula for the 547/727 race. At $2$ the two nonzero character values are $1$ and $-1$, giving $u_1(2)=0$ and $u_2(2)=2$. The weighted spline sums at centres $2$ and $4$ are therefore $0$ and $2\log(2)/3$. The prime-ideal score in the cited race omits every place above $2$, so its contributions at norms $2$ and $4$ are both zero. Under Dirichlet GRH the limiting extrema at widths $2$ and $4$ recover the two primitive traces after the explicit background correction. [Crossley, *Tuning a prime race to one prime power*](https://doi.org/10.5281/zenodo.24412652)

## Three evaluations of the spike total

1. Sum $\Lambda(n)\phi_c(\log n)/\sqrt n$ over prime powers.
2. For every $d\le\lfloor c+2h\rfloor$, sum $B((md-c)/h)$ over multiples of $d$ and give that sum coefficient $-\mu(d)\log(d)/\sqrt c$.
3. Evaluate the same lattice sums in Fourier space. With $\mathrm{sinc}(a)=\sin(a)/a$ and $\mathrm{sinc}(0)=1$, the transform of $B$ is $\mathrm{sinc}^4(\pi\xi)$ under the convention $e^{-2\pi i\xi x}$. Poisson summation gives

$$
\sum_{m\in\mathbb Z}B((md-c)/h)
=\frac hd\left[1+2\sum_{j\ge1}\mathrm{sinc}^4(\pi h j/d)\cos(2\pi c j/d)\right].
$$

The code keeps $1\le j\le K$ along with their negative-frequency partners and the zero mode. Using $|\sin a|\le1$ and $\sum_{j>K}j^{-4}\le1/(3K^3)$ bounds the entire omitted Fourier tail, after divisor weighting, by

```math
\frac{2}{3\pi^4 h^3 K^3\sqrt c}
\sum_{d\le c+2h}|\mu(d)|\,(\log d)\,d^3.
```

The CLI reports this as `fourier_tail_bound`. It bounds the Fourier truncation error at the selected center in exact arithmetic. The sum of absolute values sacrifices cancellation between divisor layers. [Poisson summation and periodization](https://www2.math.upenn.edu/~gressman/analysis/14-poisson.html)

## Subtracting the background

Write $v=(x-c)/h$ and substitute $x=e^u$ into $`\int 2\cosh(u/2)\phi_c(u)\,du`$ over the window. This gives

```math
b_c=\frac{h}{\sqrt c}\left(1+\frac1c\int_{-2}^2\frac{B(v)}{1+rv}\,dv\right).
```

Evenness turns the integral into $`\int B(v)/(1-r^2v^2)\,dv`$. The second moment $`\int v^2 B(v)\,dv=1/3`$ is the sum of four variances $1/12$. Since $|v|\le2$, the integral lies between $1+r^2/3$ and $1+r^2/[3(1-4r^2)]$. This yields the reported `background_interval`. Subtracting its endpoints in reverse order from the direct prime total yields `smoothed_discrepancy_interval`. The interval formulas are analytic bounds and the printed endpoints are rounded in floating point.

At $c=64$, $r=1/8$, and $K=64$, direct enumeration gives **1.040957942** for the spike total. The Fourier total is **1.040961696**, a gap of about $3.75\times10^{-6}$, within the tail bound of $1.71\times10^{-4}$. The smoothed discrepancy is about **0.02525**. These numbers concern this translated window. A bound for the sharp discrepancy $D(U)$ would also need control as the cutoff is recovered.
