"""A de Bruijn assembler for the abundant fraction of an RNA-Seq run.

Why this exists: every sweep in this repo needed contigs, and getting them
meant MEGAHIT or rnaSPAdes, neither of which runs on Windows. That made the
pipeline unrunnable end to end on the machine it was being developed on,
which is how Phase 3 ended up using whichever dataset happened to be
reachable rather than the right one.

This is deliberately not a general-purpose metagenome assembler. It targets
the one case this pipeline cares about: small (roughly 200-2000 nt),
high-copy-number RNA elements. Viroids and obelisks replicate, so they are
present at coverage far above ordinary host transcripts. On SRR11060618,
Peach latent mosaic viroid sits at 17,000x while host contigs run 2x to
595x. That gap is the design point.

The consequence is a simplification that would be wrong in a real assembler
and is right here: k-mers below `min_count` are discarded before the graph
is built. That throws away the entire low-abundance transcriptome, which is
what makes the memory footprint tractable in pure Python, and it keeps
exactly the fraction a replicating element lives in. Do not use this to
assemble a transcriptome. Use it to pull the abundant replicating fraction
out of one.

Assembly is strand-specific: k-mers are not canonicalised. That is
deliberate. The signature this pipeline hunts includes sense/antisense
co-occurrence, so collapsing the two strands into one node would destroy the
evidence the sweep is looking for.

Known limitation, measured rather than suspected. On a Streptococcus
sanguinis SK36 transcriptome this assembler reaches an N50 near 200 nt at
every k and min_count tried, and it fails to recover Obelisk-S.s, a 1137 nt
circular RNA documented as highly abundant in that strain. It does recover
Peach latent mosaic viroid, 337 nt at 17,349x, which is close to the easiest
possible case: short enough for a greedy walk to close and abundant enough
to dominate every branch point.

Treat this as usable for short, high-copy elements and unreliable much past
a few hundred nucleotides. Where rnaSPAdes or MEGAHIT will run, use those
and feed the FASTA to the sweep instead. See
results/phase6_obelisk_control_report.md.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from .seqio import Record, revcomp

_BASES = "ACGT"


@dataclass
class AssemblyStats:
    n_reads: int = 0
    n_kmers_total: int = 0
    n_kmers_kept: int = 0
    k: int = 31
    min_count: int = 0
    n_contigs: int = 0
    total_bp: int = 0
    n50: int = 0
    max_len: int = 0
    n_circular: int = 0

    def summary(self) -> str:
        return (
            f"reads={self.n_reads} k={self.k} min_count={self.min_count} "
            f"kmers {self.n_kmers_total}->{self.n_kmers_kept} kept | "
            f"contigs={self.n_contigs} total={self.total_bp}bp "
            f"N50={self.n50} max={self.max_len} circular={self.n_circular}"
        )


@dataclass
class Assembly:
    contigs: list[Record] = field(default_factory=list)
    stats: AssemblyStats = field(default_factory=AssemblyStats)


def count_kmers(
    reads,
    k: int = 31,
    min_count: int = 1,
    max_table: int | None = None,
) -> tuple[dict[str, int], int, int]:
    """Count k-mers across reads, then drop those below min_count.

    Reads containing anything other than ACGT contribute only the k-mers
    that avoid the ambiguous positions, so an N does not silently become a
    base.

    max_table caps peak memory. A real RNA-Seq subset produces tens of
    millions of distinct k-mers, nearly all of them sequencing errors seen
    once, and a Python dict of that size will not fit in RAM. When the table
    exceeds max_table, singletons are dropped and counting continues.

    That is lossy in a specific and acceptable direction: a k-mer is only
    discarded if it has been seen exactly once so far, so the elements this
    assembler is for, which are present at tens to thousands of copies, are
    effectively never lost, while errors are. A low-abundance transcript
    can be, which is why this must not be used as a general transcriptome
    assembler. The pruning floor rises only if one pass is not enough.
    """
    counts: Counter[str] = Counter()
    n_reads = 0
    pruned = 0
    floor = 1
    valid = set(_BASES)
    for seq in reads:
        seq = seq.strip().upper()
        n_reads += 1
        if len(seq) < k:
            continue
        # Validate the read once rather than every k-mer. Checking each k-mer
        # individually costs four string scans per position, which dominated
        # runtime and put deep subsets out of reach; a clean read means every
        # one of its k-mers is clean.
        if valid.issuperset(seq):
            segments = (seq,)
        else:
            segments = tuple(
                s for s in "".join(c if c in valid else " " for c in seq).split() if len(s) >= k
            )
        for segment in segments:
            counts.update(segment[i : i + k] for i in range(len(segment) - k + 1))
        if max_table is not None and len(counts) > max_table:
            before = len(counts)
            counts = Counter({km: c for km, c in counts.items() if c > floor})
            pruned += before - len(counts)
            # if a pass barely helped, the data is deeper than the cap and the
            # floor has to rise or we will prune on every read from here on
            if len(counts) > max_table * 0.9:
                floor += 1
    total = len(counts) + pruned
    if min_count > 1:
        counts = Counter({kmer: c for kmer, c in counts.items() if c >= min_count})
    return dict(counts), n_reads, total


def _successors(kmer: str, graph: dict[str, int]) -> list[str]:
    suffix = kmer[1:]
    return [nxt for b in _BASES if (nxt := suffix + b) in graph]


def _predecessors(kmer: str, graph: dict[str, int]) -> list[str]:
    prefix = kmer[:-1]
    return [prv for b in _BASES if (prv := b + prefix) in graph]


def build_unitigs(graph: dict[str, int], k: int) -> list[str]:
    """Walk maximal non-branching paths through the de Bruijn graph.

    A path is extended while the current node has exactly one successor and
    that successor has exactly one predecessor. Anything else is a branch or
    a junction, and stopping there is what keeps distinct molecules from
    being fused into a chimera.

    Cycles are handled explicitly: a circular molecule has no branch to stop
    at, so extension is capped by a visited set and the resulting contig
    carries the wrap-around, which is exactly the terminal repeat
    `circularity.find_circularity` looks for.
    """
    visited: set[str] = set()
    unitigs: list[str] = []

    def extend_forward(start: str) -> list[str]:
        path = [start]
        seen_local = {start}
        node = start
        while True:
            succ = _successors(node, graph)
            if len(succ) != 1:
                break
            nxt = succ[0]
            if len(_predecessors(nxt, graph)) != 1:
                break
            if nxt in seen_local or nxt in visited:
                # closed a cycle; include the wrap so the repeat is visible
                path.append(nxt)
                break
            path.append(nxt)
            seen_local.add(nxt)
            node = nxt
        return path

    # Start from nodes that cannot be extended backwards, so linear stretches
    # are emitted whole; anything left over afterwards is a pure cycle.
    starts = [km for km in graph if len(_predecessors(km, graph)) != 1
              or len(_successors(_predecessors(km, graph)[0], graph)) != 1]

    for start in starts:
        if start in visited:
            continue
        path = extend_forward(start)
        for node in path:
            visited.add(node)
        unitigs.append(path[0] + "".join(node[-1] for node in path[1:]))

    for km in graph:
        if km in visited:
            continue
        path = extend_forward(km)
        for node in path:
            visited.add(node)
        unitigs.append(path[0] + "".join(node[-1] for node in path[1:]))

    return unitigs


def remove_tips(graph: dict[str, int], k: int, max_tip: int | None = None) -> int:
    """Delete short dead-end branches. Returns the number of k-mers removed.

    A sequencing error near the end of a read creates a path that leaves the
    backbone and stops after a few nodes. Those tips are what make an
    otherwise clean graph branch, and every branch halts unitig extension.
    """
    max_tip = max_tip or 2 * k
    removed = 0
    changed = True
    while changed:
        changed = False
        for kmer in list(graph):
            if kmer not in graph:
                continue
            if _successors(kmer, graph) or len(_predecessors(kmer, graph)) != 1:
                continue
            # dead end: walk back while the path is linear
            path = [kmer]
            node = kmer
            while len(path) <= max_tip:
                preds = _predecessors(node, graph)
                if len(preds) != 1:
                    break
                prev = preds[0]
                if len(_successors(prev, graph)) != 1:
                    break  # prev is the branch point the tip hangs off
                path.append(prev)
                node = prev
            if len(path) <= max_tip:
                for node in path:
                    graph.pop(node, None)
                    removed += 1
                changed = True
    return removed


def greedy_contigs(
    graph: dict[str, int],
    k: int,
    min_ratio: float = 0.05,
) -> list[str]:
    """Trace the highest-coverage path through the graph from each seed.

    Bubble popping assumes bubbles are isolated and reconverge quickly. In a
    replicating RNA population they are neither: variants sit at a large
    fraction of positions, so the bubbles nest and a popper that only
    handles the simple case leaves the backbone in pieces. On the
    SRR11060618 subset, strict extension shredded a 337 nt viroid into
    unitigs of 25 to 73 nt whether or not bubbles were popped first.

    Greedy extension sidesteps the whole problem. The consensus genome is by
    definition the most abundant sequence at every position, so following
    the deepest successor at each step traces it directly, however tangled
    the minority variants are around it.

    Seeds are taken in descending coverage order, so the most abundant
    molecule is reconstructed first and cannot be truncated by an earlier,
    shallower path claiming its nodes. Extension stops when the next node
    falls below min_ratio of the current path's running mean, which is what
    keeps two different molecules that happen to share a k-mer from being
    fused into a chimera.
    """
    visited: set[str] = set()
    contigs: list[str] = []

    def best_next(node: str, mean: float, forward: bool) -> str | None:
        options = _successors(node, graph) if forward else _predecessors(node, graph)
        if not options:
            return None
        nxt = max(options, key=lambda o: graph[o])
        if graph[nxt] < mean * min_ratio:
            return None
        return nxt

    for seed in sorted(graph, key=lambda km: -graph[km]):
        if seed in visited:
            continue
        path = [seed]
        seen = {seed}
        total = graph[seed]

        # Crossing a stretch another contig already claimed is legitimate:
        # molecules share short repeats, and a quasispecies shares far more.
        # Re-tracing an entire existing contig is not, and blocking claimed
        # nodes outright truncated any molecule that crossed one. Obelisk-S.s
        # is 1-13% of reads in S. sanguinis SK36 and came out as fragments of
        # 194, 282 and 566 nt at 10,000-58,000x instead of one 1137 nt circle;
        # every fragment still had successors in the graph, so extension had
        # stopped with somewhere to go.
        #
        # So allow crossing, but only for a bounded run of claimed nodes.
        # That keeps the walk linear in graph size instead of letting every
        # seed re-traverse everything.
        max_revisit = 2 * k
        run = 0

        node = seed
        while True:
            nxt = best_next(node, total / len(path), forward=True)
            if nxt is None:
                break
            if nxt in seen:
                path.append(nxt)  # closed a cycle: keep the wrap as a terminal repeat
                break
            run = run + 1 if nxt in visited else 0
            if run > max_revisit:
                break
            path.append(nxt)
            seen.add(nxt)
            total += graph[nxt]
            node = nxt

        run = 0
        node = seed
        while True:
            prv = best_next(node, total / len(path), forward=False)
            if prv is None or prv in seen:
                break
            run = run + 1 if prv in visited else 0
            if run > max_revisit:
                break
            path.insert(0, prv)
            seen.add(prv)
            total += graph[prv]
            node = prv

        visited.update(path)
        contigs.append(path[0] + "".join(n[-1] for n in path[1:]))

    # Extension no longer avoids claimed nodes, so two seeds on the same
    # molecule can emit the same sequence or one contained in the other.
    # Keep the longest of each such group.
    contigs.sort(key=len, reverse=True)
    kept: list[str] = []
    for contig in contigs:
        if not any(contig in longer for longer in kept):
            kept.append(contig)
    return kept


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2
    run = 0
    for length in ordered:
        run += length
        if run >= half:
            return length
    return ordered[-1]


def assemble(
    reads,
    k: int = 31,
    min_count: int = 5,
    min_contig_len: int = 200,
    max_table: int | None = 20_000_000,
    clean: bool = True,
    greedy: bool = True,
) -> Assembly:
    """Assemble the abundant fraction of a read set into contigs.

    clean removes error tips before extension. greedy traces the deepest
    path through the graph rather than only unambiguous stretches; leave it
    on for anything that varies within a sample, which includes every
    replicating RNA element this pipeline is aimed at. Turn it off to get
    strict unitigs, which are more conservative and much more fragmented.
    """
    graph, n_reads, n_total = count_kmers(reads, k=k, min_count=min_count, max_table=max_table)
    stats = AssemblyStats(
        n_reads=n_reads, n_kmers_total=n_total, n_kmers_kept=len(graph),
        k=k, min_count=min_count,
    )
    if not graph:
        return Assembly(contigs=[], stats=stats)

    if clean:
        remove_tips(graph, k)
        stats.n_kmers_kept = len(graph)

    raw = greedy_contigs(graph, k) if greedy else build_unitigs(graph, k)
    unitigs = [u for u in raw if len(u) >= min_contig_len]
    unitigs.sort(key=len, reverse=True)

    contigs: list[Record] = []
    for i, seq in enumerate(unitigs):
        # mean k-mer coverage, the closest analogue to the assembler-reported
        # multiplicity the rest of the pipeline reads out of contig headers
        kms = [seq[j : j + k] for j in range(len(seq) - k + 1)]
        cov = sum(graph.get(km, 0) for km in kms) / max(len(kms), 1)
        contigs.append(Record(f"contig_{i} flag=1 multi={cov:.4f} len={len(seq)}", seq))

    lengths = [len(c.seq) for c in contigs]
    stats.n_contigs = len(contigs)
    stats.total_bp = sum(lengths)
    stats.n50 = _n50(lengths)
    stats.max_len = max(lengths) if lengths else 0
    return Assembly(contigs=contigs, stats=stats)


def iter_fastq(path: str, limit: int | None = None):
    """Yield sequences from a FASTQ, transparently handling .gz."""
    import gzip
    import io

    opener = (lambda p: io.TextIOWrapper(gzip.open(p, "rb"), encoding="utf-8", errors="replace")) \
        if path.endswith(".gz") else (lambda p: open(p, encoding="utf-8", errors="replace"))
    with opener(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.strip()
                if limit is not None and (i // 4) + 1 >= limit:
                    return
