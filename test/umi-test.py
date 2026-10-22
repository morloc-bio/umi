# Generated reads and UMI sets with an independent oracle for umi.
import os
import random
import tempfile


def write_pairs(n, seed):
    rng = random.Random(seed)
    d = tempfile.mkdtemp(prefix="umi-test-")
    p = os.path.join(d, "r")
    with open(p + "_1.fq", "w") as f1, open(p + "_2.fq", "w") as f2:
        for i in range(n):
            s1 = "".join(rng.choice("ACGT") for _ in range(30))
            s2 = "".join(rng.choice("ACGT") for _ in range(25))
            q1 = "".join(chr(33 + rng.randint(2, 40)) for _ in s1)
            q2 = "".join(chr(33 + rng.randint(2, 40)) for _ in s2)
            f1.write("@r%d 1:N:0\n%s\n+\n%s\n" % (i, s1, q1))
            f2.write("@r%d 2:N:0\n%s\n+\n%s\n" % (i, s2, q2))
    return p


def _read(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def check_extract(prefix, pattern, out1, out2):
    """'' when OUT holds the inputs with the UMI moved as umi_tools does."""
    for (h1, s1, q1), (h2, s2, q2), (g1, t1, r1), (g2, t2, r2) in zip(
        _read(prefix + "_1.fq"), _read(prefix + "_2.fq"), _read(out1), _read(out2)
    ):
        umi = "".join(c for c, p in zip(s1, pattern) if p == "N")
        keep = "".join(c for c, p in zip(s1, pattern) if p == "X") + s1[len(pattern):]
        keepq = "".join(c for c, p in zip(q1, pattern) if p == "X") + q1[len(pattern):]
        w1 = h1.split(" ", 1)
        w2 = h2.split(" ", 1)
        want1 = w1[0] + "_" + umi + (" " + w1[1] if len(w1) > 1 else "")
        want2 = w2[0] + "_" + umi + (" " + w2[1] if len(w2) > 1 else "")
        if (g1, t1, r1, g2, t2, r2) != (want1, keep, keepq, want2, s2, q2):
            return "read %s: got %r" % (h1, (g1, t1, g2))
    return ""


def umi_sets(seed, n):
    """n UMI sets: true UMIs with errors one base away at lower counts."""
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        counts = {}
        for _ in range(rng.randint(1, 5)):
            u = "".join(rng.choice("ACGT") for _ in range(6))
            counts[u] = counts.get(u, 0) + rng.randint(1, 50)
            for _ in range(rng.randint(0, 4)):
                k = rng.randrange(6)
                e = u[:k] + rng.choice([c for c in "ACGT" if c != u[k]]) + u[k + 1:]
                counts[e] = counts.get(e, 0) + rng.randint(1, 5)
        out.append(sorted(counts.items()))
    return out


def directional(items):
    counts = dict(items)
    umis = sorted(counts, key=lambda u: (-counts[u], u))
    adj = {u: [v for v in umis if v != u and sum(a != b for a, b in zip(u, v)) == 1 and counts[u] >= 2 * counts[v] - 1] for u in umis}
    found, comps = set(), []
    for u in umis:
        if u in found:
            continue
        comp, todo = {u}, [u]
        while todo:
            x = todo.pop()
            for y in adj[x]:
                if y not in comp:
                    comp.add(y)
                    todo.append(y)
        found |= comp
        comps.append(comp)
    return comps


def check_clusters(sets, got):
    for items, g in zip(sets, got):
        want = sorted(sorted(c) for c in directional(items))
        if sorted(sorted(c) for c in g) != want:
            return "%r: got %r, want %r" % (items, g, want)
    return ""


def fragments(seed, n):
    """(keys, UMIs, stream indices) of n fragments: 20 positions, UMIs with
    errors, some fragments without a key; keys as dup's DupKeys, the
    indices a permutation of 0 .. n-1."""
    rng = random.Random(seed)
    base = ["".join(rng.choice("ACGT") for _ in range(6)) for _ in range(8)]
    keyed, words, umis = [], [], []
    for _ in range(n):
        k = None if rng.random() < 0.05 else rng.randrange(20)
        keyed.append(k is not None)
        if k is not None:
            words += [(1 << 62) | k, 7 * k]
        u = rng.choice(base)
        if rng.random() < 0.1:
            j = rng.randrange(6)
            u = u[:j] + rng.choice([c for c in "ACGT" if c != u[j]]) + u[j + 1:]
        umis.append(u)
    ords = list(range(n))
    rng.shuffle(ords)
    return ({"keyed": keyed, "words": words}, umis, ords)


def _groups(frags):
    keys, umis, ords = frags
    groups, w = {}, 0
    for i, k in enumerate(keys["keyed"]):
        if k:
            key = (keys["words"][w], keys["words"][w + 1])
            w += 2
            groups.setdefault(key, []).append((umis[i], ords[i]))
    return groups


def dedup(frags):
    groups = _groups(frags)
    total = sum(len(g) for g in groups.values())
    mol = 0
    for g in groups.values():
        counts = {}
        for u, _ in g:
            counts[u] = counts.get(u, 0) + 1
        mol += len(directional(list(counts.items())))
    return (total, mol)


def kept(frags):
    """Per fragment, whether it is kept: a fragment without a key always is
    (nothing to deduplicate it by), a keyed one when it is its molecule's
    first fragment in the stream."""
    keep = set()
    for g in _groups(frags).values():
        counts, first = {}, {}
        for u, o in g:
            counts[u] = counts.get(u, 0) + 1
            first[u] = min(first.get(u, o), o)
        found = set()
        for comp in directional(list(counts.items())):
            own = comp - found
            found |= comp
            keep.add(min(first[u] for u in own))
    return [(not k) or o in keep for k, o in zip(frags[0]["keyed"], frags[2])]


def tmp_dir():
    import tempfile
    return tempfile.mkdtemp(prefix="umi-test-")


def same_key(n):
    return {"keyed": [True] * n, "words": [1 << 62, 5] * n}
