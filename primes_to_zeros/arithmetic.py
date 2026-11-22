"""The positive half of the logarithmic prime-power measure."""

from itertools import compress
import math


def primes_through(limit: int):
    """Yield primes up to the inclusive integer cutoff."""
    if limit < 2:
        return
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\x00\x00"
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start::p] = bytes((limit - start) // p + 1)
    yield from compress(range(limit + 1), flags)


def prime_powers(limit: int):
    """Yield (p^k, k, log p), grouped by p. Integer powers decide boundary cases."""
    if type(limit) is not int or limit < 1:
        raise ValueError("cutoff must be a positive integer")
    for p in primes_through(limit):
        n, k = p, 1
        log_p = math.log(p)
        while n <= limit:
            yield n, k, log_p
            n *= p
            k += 1


def split_powers(limit: int) -> dict:
    """Split S into first powers, squares and all exponents k >= 3."""
    weights = [[], [], []]
    for n, k, log_p in prime_powers(limit):
        bucket = min(k, 3) - 1
        weights[bucket].append(log_p / math.sqrt(n))
    parts = [math.fsum(group) for group in weights]
    spikes = math.fsum(parts)
    background = 2 * ((limit - 1) / math.sqrt(limit))
    return {
        "cutoff": limit,
        "groups": {
            name: {"count": len(group), "weight": value}
            for name, group, value in zip(("primes", "squares", "higher"), weights, parts)
        },
        "spike_total": spikes,
        "background": background,
        "discrepancy": spikes - background,
        "prime_only_discrepancy": parts[0] - background,
    }


def divisor_layers(limit: int):
    """Yield (d, μ(d), -μ(d) log d) for each nonzero divisor layer."""
    mu = [1] * (limit + 1)
    for p in primes_through(limit):
        for n in range(p, limit + 1, p):
            mu[n] = -mu[n]
        for n in range(p * p, limit + 1, p * p):
            mu[n] = 0
    for d in range(2, limit + 1):
        if mu[d]:
            yield d, mu[d], -mu[d] * math.log(d)


def at_integer(n: int) -> dict:
    """Reconstruct the log-weight at one integer and its normalized spike."""
    if type(n) is not int or n < 1:
        raise ValueError("n must be a positive integer")
    terms = [{"d": d, "mu": mu, "coefficient": weight}
             for d, mu, weight in divisor_layers(n) if n % d == 0]
    # Combine the logarithms as an integer ratio before evaluating their sum.
    numerator = math.prod(term["d"] for term in terms if term["mu"] == -1)
    denominator = math.prod(term["d"] for term in terms if term["mu"] == 1)
    base, remainder = divmod(numerator, denominator)
    if remainder:
        raise ArithmeticError("Möbius product should be 1 or a prime")
    total = math.log(base)
    return {"n": n, "log_n": math.log(n), "layers": terms, "log_argument": base, "mangoldt": total,
            "spike_weight": total / math.sqrt(n)}
