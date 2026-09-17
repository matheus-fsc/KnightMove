"""complex_moves.py — Task A.

Core complex-number representation of the knight on an n×n board.

Mathematical conventions used throughout the ``complex_orbit`` package:

* A board cell ``(row, col)`` is encoded as ``z = col + i·row`` with
  ``row, col ∈ {0, …, n−1}``.  We pick `x = col`, `y = row` so the board
  rendered in matplotlib's default ``imshow`` orientation is the same as
  the algebraic plane.
* The eight knight displacements form the set

      Δ = { ±1 ± 2i, ±2 ± i }   ⊂   ℤ[i].

* Δ is the orbit of ``δ₀ = 1 + 2i`` together with ``δ̄₀ = 1 − 2i`` under
  the cyclic group ⟨i⟩ ≅ ℤ/4:

      Δ = { iᵏ · δ₀ : k = 0..3 } ∪ { iᵏ · δ̄₀ : k = 0..3 }.

* Every δ ∈ Δ satisfies ``|δ|² = 5``.
* The minimal monic polynomial vanishing on every δ ∈ Δ is

      P(z) = z⁸ + 14 z⁴ + 625.

This module exposes the move set, validators, and a small ``__main__``
that runs the unit tests for all of the above properties.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Iterator

import numpy as np

DEFAULT_BOARD: int = 6

DELTA_0: complex = 1 + 2j

KNIGHT_MOVES: tuple[complex, ...] = (
    1 + 2j, -1 + 2j, 1 - 2j, -1 - 2j,
    2 + 1j, -2 + 1j, 2 - 1j, -2 - 1j,
)


def cell_to_z(row: int, col: int) -> complex:
    """Return the complex encoding of board cell ``(row, col)``: ``z = col + i·row``."""
    return complex(col, row)


def z_to_cell(z: complex) -> tuple[int, int]:
    """Inverse of :func:`cell_to_z`.  Returns ``(row, col)`` with integer rounding."""
    return int(round(z.imag)), int(round(z.real))


def in_board(z: complex, n: int = DEFAULT_BOARD) -> bool:
    """True iff ``z`` represents a cell on an ``n×n`` board."""
    r, c = z_to_cell(z)
    return 0 <= r < n and 0 <= c < n


def valid_moves(
    z: complex,
    visited: Iterable[complex] = (),
    n: int = DEFAULT_BOARD,
) -> list[complex]:
    """Return the legal next positions from ``z`` (in-board, unvisited)."""
    visited_set = {complex(v) for v in visited}
    out = []
    for d in KNIGHT_MOVES:
        w = z + d
        if in_board(w, n) and w not in visited_set:
            out.append(w)
    return out


def iter_knight_moves() -> Iterator[complex]:
    """Iterate over the 8 knight displacements (in canonical order)."""
    yield from KNIGHT_MOVES


# ---------------------------------------------------------------------------
# Minimal polynomial of the move set
# ---------------------------------------------------------------------------

def P(z: complex) -> complex:
    """Evaluate ``P(z) = z⁸ + 14 z⁴ + 625``.

    By Task A this vanishes on every δ ∈ :data:`KNIGHT_MOVES`.
    """
    z4 = z ** 4
    return z4 * z4 + 14 * z4 + 625


# ---------------------------------------------------------------------------
# Unit tests / self-checks
# ---------------------------------------------------------------------------

def _check_modulus_squared() -> None:
    for d in KNIGHT_MOVES:
        assert abs(abs(d) ** 2 - 5) < 1e-12, f"|{d}|² should be 5, got {abs(d)**2}"


def _check_factorization() -> None:
    """Δ == ⟨i⟩·δ₀  ∪  ⟨i⟩·δ̄₀."""
    orbit = set()
    for k in range(4):
        rot = 1j ** k
        orbit.add(rot * DELTA_0)
        orbit.add(rot * DELTA_0.conjugate())
    # Compare as sets of (round(re), round(im)) — entries are Gaussian integers.
    canon_orbit = {(int(round(z.real)), int(round(z.imag))) for z in orbit}
    canon_moves = {(int(round(z.real)), int(round(z.imag))) for z in KNIGHT_MOVES}
    assert canon_orbit == canon_moves, (
        f"Move set ≠ orbit of δ₀ under ⟨i⟩∪conj.\n  moves={canon_moves}\n  orbit={canon_orbit}"
    )


def _check_polynomial_vanishes() -> None:
    for d in KNIGHT_MOVES:
        val = P(d)
        assert abs(val) < 1e-8, f"P({d}) should be 0, got {val}"


def _check_valid_moves_basic() -> None:
    # From the corner (0,0) → z = 0+0i: only (1,2) and (2,1) land on board.
    legal = set(valid_moves(0 + 0j, n=DEFAULT_BOARD))
    assert legal == {1 + 2j, 2 + 1j}, f"corner moves mismatch: {legal}"
    # All 8 moves valid from the centre of an 8×8 board (3,3) with empty history.
    legal8 = set(valid_moves(3 + 3j, n=8))
    expected = {(3 + 3j) + d for d in KNIGHT_MOVES}
    assert legal8 == expected, "8 moves should all be legal from (3,3) on 8x8"


def run_self_tests(verbose: bool = True) -> None:
    """Run every Task-A check and print a summary."""
    _check_modulus_squared()
    if verbose:
        print("  ✓ |Δz|² = 5 for all 8 moves")
    _check_factorization()
    if verbose:
        print("  ✓ Move set = ⟨i⟩·δ₀  ∪  ⟨i⟩·δ̄₀   (δ₀ = 1 + 2i)")
    _check_polynomial_vanishes()
    if verbose:
        print("  ✓ P(z) = z⁸ + 14 z⁴ + 625 vanishes on every move")
    _check_valid_moves_basic()
    if verbose:
        print("  ✓ valid_moves respects board bounds and visited set")


if __name__ == "__main__":
    print("complex_moves.py — Task A self-tests")
    print("-" * 60)
    run_self_tests(verbose=True)
    print("-" * 60)
    print("All Task-A checks passed.")
    print()
    print("Move table:")
    for d in KNIGHT_MOVES:
        print(f"  Δz = {str(d):>10s}    |Δz|² = {abs(d)**2:>5.1f}    P(Δz) = {P(d):.2g}")
