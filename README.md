# Packing the traces of rectangles on horizontal lines

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

Let F be a finite family of closed axis-parallel rectangles in the plane and L a finite set of horizontal
lines such that every rectangle meets some line. Intersecting each rectangle with the union of the lines gives
a family of p-intervals. McGinnis and Zerbib conjectured that the packing number of the traces is at most a
constant times the packing number of F. The paper proves the bound 3 nu(F) - 1 for every nonempty family and
every such set of lines, and (3 - epsilon) nu(F) for an absolute epsilon > 0. For a prescribed set of lines
the constant cannot be below 5/2, and 5/2 is the best constant for four or five lines.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper (build with `pdflatex main.tex` twice).
- `code/check_traces.py`: checks the lemmas of Sections 2 and 3 on random families of rectangles with small
  integer coordinates (three generators: free placement, stacked on consecutive line blocks, and dense stacks
  aimed at four-cycles), then the five-cycle example and its variants, and rebuilds the Thomas-Walls graph
  T'_3. Run `python code/check_traces.py SEED TRIALS`; it needs `networkx`.
- `code/check_masks.py`: the pair coverage of the lists of sets in Section 4.

With seeds 4 and 5 and 20,000 trials each, `check_traces.py` reported no failure over 72,625 graph edges,
439 four-cycles and 109,680 pairs of representatives with the same block.

`SHA256SUMS`: sha256 of each file above.
