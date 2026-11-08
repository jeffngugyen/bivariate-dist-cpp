"""Compare the real CLI/Python API with an independent exhaustive oracle.

The oracle enumerates bit strings directly; it does not use Gaussian
elimination or the solver's row-reduction algorithm.
"""
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile
import unittest

BINARY = Path(os.environ.get("DIST_BINARY", "build/dist")).resolve()
BACKEND = os.environ.get("DIST_BACKEND", "cli")
if BACKEND == "python":
    import dist_py


def span(rows):
    result = {0}
    for row in rows:
        result |= {vector ^ row for vector in tuple(result)}
    return result


def kernel(rows, width):
    return {vector for vector in range(1 << width)
            if all((vector & row).bit_count() % 2 == 0 for row in rows)}


def oracle(x, z, width):
    first = kernel(x, width) - span(z)
    second = kernel(z, width) - span(x)
    return min(v.bit_count() for v in first), min(v.bit_count() for v in second)


def rows_from_masks(masks, width):
    return [[(mask >> bit) & 1 for bit in range(width)] for mask in masks]


class SolverTests(unittest.TestCase):
    def cli(self, x, z, l, m):
        with tempfile.TemporaryDirectory() as directory:
            paths = [Path(directory) / name for name in ("X.txt", "Z.txt")]
            for path, rows in zip(paths, (x, z)):
                path.write_text("".join(" ".join(map(str, row)) + "\n" for row in rows))
            return subprocess.run([str(BINARY), *map(str, paths), str(l), str(m)],
                                  text=True, capture_output=True, timeout=10)

    def solve(self, x, z, l, m):
        if BACKEND == "python":
            return dist_py.compute_dist(x, z, l, m)
        run = self.cli(x, z, l, m)
        self.assertEqual(run.returncode, 0, run.stderr)
        match = re.fullmatch(r"d_Z=(\d+) d_X=(\d+) distance=(\d+)\n", run.stdout)
        self.assertIsNotNone(match, run.stdout)
        first, second, minimum = map(int, match.groups())
        self.assertEqual(minimum, min(first, second))
        return first, second

    def rejected(self, x, z, l, m, message):
        if BACKEND == "python":
            with self.assertRaisesRegex(ValueError, message):
                dist_py.compute_dist(x, z, l, m)
        else:
            run = self.cli(x, z, l, m)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn(message, run.stderr)

    def test_four_qubit_known_answer(self):
        self.assertEqual(self.solve([[1, 1, 1, 1]], [[1, 1, 1, 1]], 1, 2), (2, 2))

    def test_asymmetric_distances_preserve_tuple_order(self):
        # Z checks: Z1Z2, Z2Z3, Z3Z4. Z distance 1; X distance 4.
        self.assertEqual(self.solve([], rows_from_masks([3, 6, 12], 4), 1, 2), (1, 4))
        self.assertEqual(self.solve(rows_from_masks([3, 6, 12], 4), [], 1, 2), (4, 1))

    def test_empty_check_sets(self):
        self.assertEqual(self.solve([], [], 1, 1), (1, 1))

    def test_redundant_and_zero_rows(self):
        x = [[1, 1, 1, 1], [0, 0, 0, 0], [1, 1, 1, 1]]
        self.assertEqual(self.solve(x, x[::-1], 1, 2), (2, 2))

    def test_200_seeded_commuting_codes_against_oracle(self):
        rng = random.Random(20261108)
        for width in (2, 4, 6, 8):
            accepted = 0
            while accepted < 50:
                x = [rng.randrange(1 << width) for _ in range(rng.randrange(width + 1))]
                allowed = sorted(kernel(x, width))
                z = [rng.choice(allowed) for _ in range(rng.randrange(width + 1))]
                # rank(X)+rank(Z) == n encodes no logical qubits.
                if len(span(x)) * len(span(z)) == (1 << width):
                    continue
                expected = oracle(x, z, width)
                with self.subTest(width=width, case=accepted):
                    actual = self.solve(rows_from_masks(x, width), rows_from_masks(z, width),
                                        1, width // 2)
                    self.assertEqual(actual, expected)
                accepted += 1

    def test_noncommuting_checks(self):
        self.rejected([[1, 0]], [[1, 0]], 1, 1, "must commute")

    def test_ragged_and_wrong_width_rows(self):
        self.rejected([[1]], [], 1, 1, "columns")
        self.rejected([[1, 0], [0]], [], 1, 1, "columns")
        self.rejected([[1, 0, 0]], [], 1, 1, "columns")

    def test_nonbinary_input(self):
        self.rejected([[2, 0]], [], 1, 1, "0 or 1")
        self.rejected([[-1, 0]], [], 1, 1, "0 or 1")

    def test_bad_dimensions(self):
        for l, m in ((0, 1), (-1, 1), (1, 0), (1, 65), (2147483647, 2147483647)):
            with self.subTest(l=l, m=m):
                self.rejected([], [], l, m, "positive")

    def test_no_logical_qubits(self):
        self.rejected([[1, 0]], [[0, 1]], 1, 1, "no logical qubits")

    def test_64_dimensional_nullspace_is_rejected_before_shift(self):
        self.rejected([], [], 4, 8, "nullspace dimension < 64")
        # Also check the second nullspace before starting the first search.
        self.rejected(rows_from_masks([1], 64), [], 4, 8, "nullspace dimension < 64")

    def test_cli_missing_file_and_usage(self):
        with tempfile.TemporaryDirectory() as directory:
            missing = str(Path(directory) / "missing.txt")
            run = subprocess.run([str(BINARY), missing, missing, "1", "1"],
                                 text=True, capture_output=True, timeout=10)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("cannot open", run.stderr)
        run = subprocess.run([str(BINARY)], text=True, capture_output=True, timeout=10)
        self.assertEqual(run.returncode, 2)
        self.assertIn("Usage:", run.stderr)
        help_run = subprocess.run([str(BINARY), "--help"], text=True, capture_output=True, timeout=10)
        self.assertEqual(help_run.returncode, 0)

    def test_cli_rejects_malformed_tokens(self):
        if BACKEND != "cli":
            self.skipTest("CLI text parsing only")
        self.rejected([["hello", 0]], [], 1, 1, "0 or 1")

    @unittest.skipUnless(BACKEND == "python", "Python binding only")
    def test_python_keywords_and_defaults(self):
        self.assertEqual(dist_py.compute_dist(X=[[1, 1, 1, 1]], Z=[[1, 1, 1, 1]], l=1, m=2),
                         (2, 2))
        # The two-argument call uses l=4, m=7; malformed width fails safely.
        with self.assertRaisesRegex(ValueError, "columns"):
            dist_py.compute_dist([[1, 1]], [])


if __name__ == "__main__":
    unittest.main()
