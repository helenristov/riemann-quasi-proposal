import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

from primes_to_zeros.arithmetic import at_integer, split_powers
from primes_to_zeros.poisson import cubic_spline, fourier_lattice, moving_window


def trial_prime_power(n):
    """Independently find the first prime divisor and its exact exponent."""
    if n < 2:
        return None
    p = next((d for d in range(2, math.isqrt(n) + 1) if n % d == 0), n)
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return (p, k) if n == 1 else None


class ArithmeticTests(unittest.TestCase):
    def test_exact_cancellation_at_six_and_eight(self):
        six, eight = at_integer(6), at_integer(8)
        self.assertEqual([(t["d"], t["mu"]) for t in six["layers"]], [(2, -1), (3, -1), (6, 1)])
        self.assertEqual((six["log_argument"], six["spike_weight"]), (1, 0.0))
        self.assertEqual(eight["log_argument"], 2)
        self.assertAlmostEqual(eight["spike_weight"], math.log(2) / math.sqrt(8))

    def test_divisor_layers_against_factorization(self):
        for n in range(1, 301):
            power = trial_prime_power(n)
            self.assertEqual(at_integer(n)["log_argument"], power[0] if power else 1)

    def test_inclusive_exponent_groups(self):
        for cutoff in (1, 2, 7, 8, 9, 16, 27, 49, 242, 243, 255, 256, 400):
            weights = [[], [], []]
            for n in range(2, cutoff + 1):
                power = trial_prime_power(n)
                if power:
                    p, k = power
                    weights[min(2, k - 1)].append(math.log(p) / math.sqrt(n))
            result = split_powers(cutoff)
            for name, values in zip(("primes", "squares", "higher"), weights):
                self.assertEqual(result["groups"][name]["count"], len(values))
                self.assertAlmostEqual(result["groups"][name]["weight"], math.fsum(values))
            self.assertAlmostEqual(result["background"], 4 * math.sinh(math.log(cutoff) / 2))

    def test_million_and_reduction(self):
        result = split_powers(1_000_000)
        p, square, higher = [result["groups"][g] for g in ("primes", "squares", "higher")]
        self.assertEqual([p["count"], square["count"], higher["count"]], [78498, 168, 68])
        self.assertAlmostEqual(square["weight"], 5.609510475393043, places=12)
        self.assertAlmostEqual(higher["weight"], 2.228729918576302, places=12)
        self.assertAlmostEqual(result["discrepancy"], -3.09119118403, places=10)
        self.assertAlmostEqual(result["prime_only_discrepancy"] + square["weight"] + higher["weight"],
                               result["discrepancy"], places=10)
        self.assertLess(square["weight"], math.log(1_000_000)**2 / 4)

    def test_empty_origin_and_cutoff_domain(self):
        self.assertEqual(split_powers(1)["discrepancy"], 0.0)
        self.assertEqual(at_integer(1)["log_argument"], 1)
        for value in (0, -1, 2.5):
            with self.assertRaises(ValueError):
                split_powers(value)


class WindowTests(unittest.TestCase):
    def test_spline_mass_partition_and_knots(self):
        self.assertEqual(cubic_spline(0.0), 2/3)
        self.assertEqual(cubic_spline(-1.0), 1/6)
        self.assertEqual(cubic_spline(2.0), 0.0)
        for shift in (0.0, 0.23, 0.75, 1.0):
            self.assertAlmostEqual(sum(cubic_spline(i - shift) for i in range(-3, 5)), 1)

    def test_zero_frequency_of_aligned_lattice(self):
        for d in (2, 3, 7):
            self.assertAlmostEqual(fourier_lattice(d, 61, 2 * d, 1), 2.0)

    def test_fourier_vs_prime_and_multiples(self):
        for center, width in ((64, 0.125), (37, 0.137), (10, 0.45)):
            for modes in (7, 64):
                result = moving_window(center, width, modes)
                self.assertAlmostEqual(result["prime_side"], result["divisor_side"], places=12)
                self.assertLessEqual(abs(result["fourier_minus_divisor"]), result["fourier_tail_bound"] + 1e-13)
        reference = moving_window(64, 0.125, 64)
        self.assertAlmostEqual(reference["prime_side"], 1.0409579422467623, places=12)
        self.assertAlmostEqual(reference["fourier_side"], 1.0409616959829704, places=12)

    def test_background_interval_by_independent_quadrature(self):
        for center, r in ((3, 0.25), (64, 0.125), (101, 0.31)):
            result = moving_window(center, r, 16)
            a, b = map(math.log, result["support_x"])
            h = center * r
            # Midpoint quadrature on the original log axis, including 2*cosh(u/2).
            step = (b - a) / 6000

            def density(u):
                x = math.exp(u)
                return 2 * math.cosh(u / 2) * math.sqrt(x / center) * cubic_spline((x - center) / h)

            background = step * math.fsum(density(a + (i + 0.5) * step) for i in range(6000))
            low, high = result["background_interval"]
            self.assertLess(low, background)
            self.assertLess(background, high)
            err = result["prime_side"] - background
            lower, upper = result["smoothed_discrepancy_interval"]
            self.assertLess(lower, err)
            self.assertLess(err, upper)

    def test_fixed_shape_and_tail_scaling(self):
        a = moving_window(64, 0.125, 10)
        b = moving_window(128, 0.125, 10)
        self.assertEqual([2 * x for x in a["support_x"]], b["support_x"])
        self.assertAlmostEqual(a["fourier_tail_bound"], 8 * moving_window(64, 0.125, 20)["fourier_tail_bound"])

    def test_window_domain(self):
        for args in ((1, 0.1, 3), (5, 0, 4), (5, 0.5, 4), (5, math.nan, 4),
                     (5, 0.1, 0), (64, 1e-30, 8)):
            with self.assertRaises(ValueError):
                moving_window(*args)


class EntryPointTests(unittest.TestCase):
    def run_cli(self, *args):
        root = Path(__file__).resolve().parents[1]
        return subprocess.run([sys.executable, "-m", "primes_to_zeros", *args],
                              cwd=root, text=True, capture_output=True)

    def test_json_commands(self):
        for args, length in ((("layers", "6", "8"), 2), (("split", "256"), 1),
                             (("poisson", "64", "128", "256", "--modes", "64"), 3)):
            result = self.run_cli(*args)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(len(json.loads(result.stdout)), length)

    def test_invalid_cli_input(self):
        result = self.run_cli("poisson", "64", "--relative-width", "0.9")
        self.assertEqual(result.returncode, 2)
        self.assertIn("relative width", result.stderr)
