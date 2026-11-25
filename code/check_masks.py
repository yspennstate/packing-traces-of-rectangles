"""Pair coverage of the mask lists in the small-line table.
p lines give positions 0..p; a rectangle with line block [a,b] has the pair (a, b+1);
mask m keeps it when bit a is set and bit b+1 is not. q = least number of masks
keeping a pair i<j; the factor is len(masks)/q."""
from fractions import Fraction
rows = {3: [3, 5], 5: [7, 11, 13, 19, 21], 6: [11, 13, 19, 23, 29, 37, 39, 43]}
for p, masks in rows.items():
    q = min(sum(1 for m in masks if (m >> i) & 1 and not (m >> j) & 1)
            for i in range(p + 1) for j in range(i + 1, p + 1))
    print(p, "lines", len(masks), "masks, least coverage", q, "factor", Fraction(len(masks), q))
    # fewer lines: the same masks restricted to the first positions
    for pp in range(1, p):
        qq = min(sum(1 for m in masks if (m >> i) & 1 and not (m >> j) & 1)
                 for i in range(pp + 1) for j in range(i + 1, pp + 1))
        assert qq >= q, (p, pp, qq)
print("ok")
