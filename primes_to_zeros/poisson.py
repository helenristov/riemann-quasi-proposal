"""A translated compact window on the same measure as the sharp discrepancy."""

import math

from .arithmetic import divisor_layers, prime_powers


def cubic_spline(v: float) -> float:
    """Centered convolution of four unit intervals. Integral 1, support [-2, 2]."""
    t = abs(v)
    if t < 1:
        return (4 - 6 * t**2 + 3 * t**3) / 6
    if t < 2:
        return (2 - t)**3 / 6
    return 0.0


def fourier_lattice(d: int, center: int, width: float, modes: int) -> float:
    """Approximate the unweighted sum of spline values on dZ by Poisson summation."""
    def mode(k):
        angle = math.pi * width * k / d
        sinc = math.sin(angle) / angle
        # Reduce the rational phase before computing the cosine.
        phase = math.tau * ((center * k) % d) / d
        return sinc**4 * math.cos(phase)

    return (width / d) * (1 + 2 * math.fsum(mode(k) for k in range(1, modes + 1)))


def moving_window(center: int, relative_width: float, modes: int) -> dict:
    """Evaluate one log-translated window by prime powers, divisors, and Fourier modes."""
    if type(center) is not int or center < 1:
        raise ValueError("center must be a positive integer")
    if not (0 < relative_width < 0.5):
        raise ValueError("relative width must be strictly between 0 and 0.5")
    if type(modes) is not int or modes < 1:
        raise ValueError("modes must be a positive integer")
    width = center * relative_width
    lower, upper = center - 2 * width, center + 2 * width
    if lower < 1:
        raise ValueError("choose a window whose left edge is at least 1")
    if not lower < center < upper:
        raise ValueError("window edges must be distinguishable at floating-point precision")
    limit = math.floor(upper)
    root = math.sqrt(center)
    layers = list(divisor_layers(limit))

    # Apply the test to the marks in their original normalization.
    prime_side = math.fsum(
        (log_p / math.sqrt(n)) * math.sqrt(n / center) * cubic_spline((n - center) / width)
        for n, _, log_p in prime_powers(limit) if n >= lower
    )

    def on_multiples(d):
        lo, hi = math.ceil(lower / d), math.floor(upper / d)
        return math.fsum(cubic_spline((m * d - center) / width) for m in range(lo, hi + 1))

    divisor_side = math.fsum(weight * on_multiples(d) for d, _, weight in layers) / root
    fourier_side = math.fsum(weight * fourier_lattice(d, center, width, modes)
                             for d, _, weight in layers) / root
    # Bound both tails using |sinc(a)| <= 1/|a| and sum_{k>K} k^-4 <= 1/(3K^3).
    tail = (2 / (3 * math.pi**4 * width**3 * modes**3 * root)
            * math.fsum(abs(weight) * d**3 for d, _, weight in layers))

    # The spline is even with second moment 1/3. See docs/poisson.md.
    variance = relative_width**2 / 3
    bg_low = width / root * (1 + (1 + variance) / center)
    bg_high = width / root * (1 + (1 + variance / (1 - 4 * relative_width**2)) / center)
    return {
        "center": center, "relative_width": relative_width, "support_x": [lower, upper],
        "modes_per_side": modes,
        "prime_side": prime_side, "divisor_side": divisor_side, "fourier_side": fourier_side,
        "fourier_minus_divisor": fourier_side - divisor_side, "fourier_tail_bound": tail,
        "background_interval": [bg_low, bg_high],
        "smoothed_discrepancy_interval": [prime_side - bg_high, prime_side - bg_low],
    }
