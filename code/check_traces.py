"""Random check of the rectangle-trace claims (McGinnis-Zerbib Conjecture 3.8).

Rectangles: closed, integer coordinates on a small grid (ties and touching are common).
Lines: a random set of horizontal lines hitting every rectangle.
Checks, for every maximum disjoint trace family and its representatives:
  (1) intersecting representatives have disjoint, consecutive line blocks;
  (2) the intersection graph of the representatives is triangle-free and planar;
  (3) nu(H) <= 3 nu(F) - 1  and  nu(H) <= 4 nu(F) (the colouring bound, via the best directed cut);
  (4) no four-cycle in which both opposite pairs have different blocks.
"""
import random, itertools, sys
import networkx as nx

def inter(a, b):  # closed intervals
    return a[0] <= b[1] and b[0] <= a[1]

def rect_meet(r, s):
    return inter(r[0], s[0]) and inter(r[1], s[1])

def trace_meet(r, s, lines):
    if not inter(r[0], s[0]):
        return False
    return any(r[1][0] <= y <= r[1][1] and s[1][0] <= y <= s[1][1] for y in lines)

def mis(n, adj):
    """Maximum independent set size and one witness (exact, small n)."""
    best = [0, ()]
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    def rec(cand, chosen):
        if len(chosen) + len(cand) <= best[0]:
            return
        if not cand:
            best[0], best[1] = len(chosen), tuple(chosen)
            return
        v = cand[0]
        rec([u for u in cand[1:] if u not in adj[v]], chosen + [v])
        rec(cand[1:], chosen)
    rec(order, [])
    return best[0], best[1]

def block(r, lines):
    idx = [i for i, y in enumerate(lines) if r[1][0] <= y <= r[1][1]]
    return idx[0], idx[-1]

COUNT = {"C4": 0, "sameblock": 0, "edges": 0}

def trial(rng, n, G, stacked=False, dense=False):
    R = []
    if dense:
        # aimed at four-cycles: few lines, rectangles one or two blocks tall and narrow in x,
        # packed on a short x-range so that representatives on adjacent blocks touch often
        q = rng.randint(2, 4)
        for _ in range(n):
            a = rng.randint(0, q); b = min(q, a + rng.randint(0, 1))
            x0 = rng.randint(0, G); x1 = min(G, x0 + rng.randint(0, 2))
            R.append(((x0, x1), (2 * a - 1, 2 * b + 1)))
        ys = list(range(0, 2 * q + 1, 2))
    elif stacked:
        # lines at 0,2,..,2q; each rectangle spans a run of consecutive lines and reaches
        # exactly to the gap midpoint (odd y) so that neighbours on adjacent blocks touch
        q = rng.randint(2, 6)
        for _ in range(n):
            a = rng.randint(0, q); b = rng.randint(a, min(q, a + rng.randint(0, 2)))
            x0 = rng.randint(0, G); x1 = rng.randint(x0, min(G, x0 + rng.randint(0, G // 2)))
            R.append(((x0, x1), (2 * a - 1, 2 * b + 1)))
        ys = list(range(0, 2 * q + 1, 2))
    else:
        for _ in range(n):
            x0 = rng.randint(0, G); x1 = rng.randint(x0, min(G, x0 + rng.randint(0, G // 2)))
            y0 = rng.randint(0, G); y1 = rng.randint(y0, min(G, y0 + rng.randint(0, G // 2)))
            R.append(((x0, x1), (y0, y1)))
        ys = sorted(set(rng.sample(range(0, G + 1), rng.randint(1, G // 2 + 1))))
    # add lines until every rectangle is hit
    for r in R:
        if not any(r[1][0] <= y <= r[1][1] for y in ys):
            ys.append(rng.randint(r[1][0], r[1][1]))
    lines = sorted(set(ys))
    adjF = [set(j for j in range(n) if j != i and rect_meet(R[i], R[j])) for i in range(n)]
    adjH = [set(j for j in range(n) if j != i and trace_meet(R[i], R[j], lines)) for i in range(n)]
    nuF, _ = mis(n, adjF)
    nuH, rep = mis(n, adjH)
    bad = []
    if nuH > 3 * nuF - 1:
        bad.append("3nu-1")
    g = nx.Graph(); g.add_nodes_from(rep)
    for i, j in itertools.combinations(rep, 2):
        if rect_meet(R[i], R[j]):
            g.add_edge(i, j)
            (a1, b1), (a2, b2) = block(R[i], lines), block(R[j], lines)
            if not (b1 + 1 == a2 or b2 + 1 == a1):
                bad.append("consecutive")
    if any(len(set(g[u]) & set(g[v])) for u, v in g.edges()):
        bad.append("triangle")
    if not nx.check_planarity(g)[0]:
        bad.append("planar")
    # colouring bound: best directed cut of positions 0..p
    p = len(lines)
    best = 0
    for m in range(1 << (p + 1)):
        c = sum(1 for i in rep if (m >> block(R[i], lines)[0]) & 1 and not (m >> (block(R[i], lines)[1] + 1)) & 1)
        best = max(best, c)
        if best * 4 >= nuH:
            break
    if best * 4 < nuH:
        bad.append("cut4")
    # representatives with the same block are joined by no path with three edges
    COUNT["edges"] += g.number_of_edges()
    for u, v in itertools.combinations(rep, 2):
        if block(R[u], lines) == block(R[v], lines):
            COUNT["sameblock"] += 1
            for p1 in g[u]:
                for p2 in g[p1]:
                    if p2 not in (u, v) and v in g[p2] and p1 != v:
                        bad.append("3path")
    # four-cycles: one opposite pair has equal blocks
    for c in nx.simple_cycles(g, length_bound=4):
        if len(c) == 4:
            COUNT["C4"] += 1
            B = [block(R[v], lines) for v in c]
            if B[0] != B[2] and B[1] != B[3]:
                bad.append("C4")
    return nuF, nuH, len(lines), bad, R, lines

if __name__ == "__main__":
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
    worst = 0.0; fails = 0; hist = {}
    for t in range(T):
        n = rng.randint(3, 16); G = rng.choice([6, 8, 10, 14])
        if t % 3 == 2:
            n = rng.randint(10, 20); G = rng.choice([4, 6, 8])
        nuF, nuH, p, bad, R, lines = trial(rng, n, G, stacked=(t % 3 == 1), dense=(t % 3 == 2))
        hist[(nuF, nuH)] = hist.get((nuF, nuH), 0) + 1
        worst = max(worst, nuH / nuF)
        if bad:
            fails += 1
            if fails <= 5:
                print("FAIL", bad, R, lines)
    print("trials", T, "fails", fails, "worst ratio %.3f" % worst)
    print("checked: graph edges", COUNT["edges"], "four-cycles", COUNT["C4"], "same-block pairs", COUNT["sameblock"])
    print("ratio>=2 cases", sorted((k, v) for k, v in hist.items() if k[1] >= 2 * k[0]))
    # the five-cycle example, its doubled version, and the 1,3,5 line choice
    C5 = [((0, 4), (-1, 1)), ((0, 1), (1, 5)), ((0, 4), (5, 7)), ((3, 4), (3, 5)), ((3, 4), (1, 3))]
    two = C5 + [((x0 + 10, x1 + 10), y) for (x0, x1), y in C5]
    thick = [(x, (y0 - 0.25, y1 + 0.25)) for x, (y0, y1) in C5]
    for name, fam, L in (("C5 lines 0246", C5, [0, 2, 4, 6]), ("2xC5 lines 0246", two, [0, 2, 4, 6]),
                         ("C5 lines 135", C5, [1, 3, 5]), ("thick C5 lines 0246", thick, [0, 2, 4, 6]),
                         ("C5 lines -1/2,0,2,4,6", C5, [-0.5, 0, 2, 4, 6])):
        n = len(fam)
        aF = [set(j for j in range(n) if j != i and rect_meet(fam[i], fam[j])) for i in range(n)]
        aH = [set(j for j in range(n) if j != i and trace_meet(fam[i], fam[j], L)) for i in range(n)]
        print(name, "nuF", mis(n, aF)[0], "nuH", mis(n, aH)[0])
    # Thomas-Walls: T1 = K4 on a,b,c,d; T2: delete ab, add x,y,z; T3: delete yz, add w,u,v; T3' deletes cd and uv.
    T = nx.complete_graph("abcd")
    T.remove_edge("a", "b"); T.add_edges_from([("a", "x"), ("b", "y"), ("b", "z"), ("x", "y"), ("x", "z"), ("y", "z")])
    T.remove_edge("y", "z"); T.add_edges_from([("y", "w"), ("z", "u"), ("z", "v"), ("w", "u"), ("w", "v"), ("u", "v")])
    T.remove_edges_from([("c", "d"), ("u", "v")])
    W = nx.cycle_graph(4)  # 0-1-2-3, opposite pairs (0,2),(1,3)
    W.add_edges_from([(0, 4), (4, 5), (5, 2), (1, 6), (6, 7), (7, 3)])
    gm = nx.algorithms.isomorphism.GraphMatcher(T, W)
    print("T3' vertices", T.number_of_nodes(), "edges", T.number_of_edges(), "triangle-free", sum(nx.triangles(T).values()) == 0,
          "planar", nx.check_planarity(T)[0], "contains W", gm.subgraph_is_monomorphic())
