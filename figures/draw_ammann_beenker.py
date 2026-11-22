"""Draw the README's Ammann–Beenker patch by projecting the integer lattice Z^4.

Run this file with Python 3 to regenerate the adjacent SVG, without dependencies.
The internal window is the regular octagon of side one centered at the origin.
"""

from itertools import combinations, product
from math import ceil, sqrt
from pathlib import Path

DIAGONAL = 1 / sqrt(2)
WINDOW_APOTHEM = (1 + sqrt(2)) / 2
PATCH_APOTHEM = 5.1
SQUARE = "#218a86"
RHOMBUS = "#72508d"


def octagon_radius(z):
    """Largest projection onto the normals of an axis-aligned regular octagon."""
    x, y = abs(z.real), abs(z.imag)
    return max(x, y, (x + y) * DIAGONAL)


def step(index, axis):
    neighbor = list(index)
    neighbor[axis] += 1
    return tuple(neighbor)


def tiles(radius):
    # Physical basis directions are 0, 45, 90, 135 degrees. The internal
    # projection multiplies each angle by three, as in the octagonal star map.
    vertices = {}
    bound = ceil(radius + 2)
    for index in product(range(-bound, bound + 1), repeat=4):
        a, b, c, d = index
        u, v = (b - d) * DIAGONAL, (b + d) * DIAGONAL
        internal = complex(a - u, v - c)
        point = complex(a + u, c + v)
        if octagon_radius(internal) < WINDOW_APOTHEM and octagon_radius(point) < radius + 1:
            vertices[index] = point

    for index, p in vertices.items():
        for i, j in combinations(range(4), 2):
            along_i, along_j = step(index, i), step(index, j)
            opposite = step(along_i, j)
            if along_i in vertices and along_j in vertices and opposite in vertices:
                q = vertices[opposite]
                if octagon_radius((p + q) / 2) <= radius:
                    kind = "square" if j - i == 2 else "rhombus"
                    yield kind, (p, vertices[along_i], q, vertices[along_j])


def draw():
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="780" height="396" viewBox="0 0 780 396" role="img" aria-labelledby="title desc">',
        '<title id="title">Ammann–Beenker tiling</title>',
        '<desc id="desc">An eightfold patch of teal squares and violet 45-degree rhombi, with enlarged examples of the two tiles.</desc>',
        '<rect width="780" height="396" rx="12" fill="#0d1824"/>',
        '<g stroke="#e1eeee" stroke-width="1.15" stroke-linejoin="round">',
    ]
    for kind, corners in tiles(PATCH_APOTHEM):
        points = " ".join(f"{221 + 32*z.real:.3f},{198 - 32*z.imag:.3f}" for z in corners)
        fill = SQUARE if kind == "square" else RHOMBUS
        lines.append(f'<polygon fill="{fill}" points="{points}"/>')
    lines += [
        '</g>',
        '<path d="M445 54V342" stroke="#364859"/>',
        '<g font-family="system-ui, sans-serif" fill="#edf3f5">',
        '<text x="478" y="88" font-size="25" font-weight="600">Two tile shapes</text>',
        f'<path d="M484 119h50v50h-50z" fill="{SQUARE}" stroke="#e1eeee" stroke-width="1.5"/>',
        '<text x="586" y="152" font-size="22">Square</text>',
        f'<path d="M478 259h50l35.355 -35.355h-50z" fill="{RHOMBUS}" stroke="#e1eeee" stroke-width="1.5"/>',
        '<text x="586" y="250" font-size="22">45° rhombus</text>',
        '</g>',
        '</svg>',
    ]
    Path(__file__).with_name("ammann-beenker.svg").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    draw()
