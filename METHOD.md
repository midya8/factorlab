# FactorLab method and limits

## Finite input and indexing

The input has exactly `n`, `factors`, `evidence`, `map_vars`, `order`. There are one to five binary variables and one to ten factors. A factor has exactly `scope` and `values`. Scope indices are sorted and distinct. The local table uses little-endian masks: bit j is the value of scope[j]. Values are nonnegative integer/fraction strings of at most 12 characters, with reduced numerator/denominator <=10000. Empty scope is a scalar. Duplicate scopes are allowed and multiply as distinct factors. Zero potentials are allowed; negative values are not.

Evidence maps canonical string variable indices to integer bits. MAP variables must be sorted, unique and unobserved. The order must permute every unobserved variable. Empty MAP scope means an all-sum query. All unobserved variables in MAP scope means an all-max query. All variables observed means no elimination and a product of scalars.

Native JSON rejects duplicate keys, nonfinite constants, Boolean/noninteger index fields and UTF-8 input beyond 32KiB. The browser limits input to 32,768 characters, JSON.parse keeps the last duplicate key, and integer-valued JSON numeric tokens such as `1.0` become indistinguishable from `1` and are accepted for index fields. Potentials are always rational strings. No binary64 inference is used; displayed bar widths and SVG coordinates alone use ordinary numbers.

## Elimination invariant

Evidence is applied by indexing original tables, not by renormalizing them. Each unobserved variable absent from every sliced factor receives an explicit `[1,1]` unit factor. This preserves the counting measure: summing an absent binary variable multiplies total mass by two, whereas maximizing it does not.

At each step, gather all factors involving v. Their union scope is joined by exact multiplication. Reduce the two entries for v using sum or max and replace that bucket by the resulting factor. Keep all other factors, including scalars. Nonnegative multiplication distributes over max; multiplication distributes over sum. The product of current factors is therefore the dense partial elimination of the original joint potential, even for a mixed schedule. Final scalars are multiplied, never dropped.

The independent reference enumerates every evidence-consistent full assignment directly from the original factors, not the elimination output. It computes the evidence-restricted partition Z_e, every joint weight, each MAP-assignment mass after summing hidden variables, all joint-mode ties and all marginal-MAP ties. Probabilities divide by Z_e only when Z_e>0. Z_e is not P(e) for arbitrary unnormalized potentials; P(e) would require the unconditioned normalizer as well. That extra normalizer is not computed here.

## Why mixed order changes the question

With nonnegative weights, all sums commute, all maxima commute, but sum and max do not generally commute. For any finite table g, `max_a sum_b g(a,b) <= sum_b max_a g(a,b)`. Each column maximum dominates every fixed-row entry. Sum the inequality and then maximize the row. Moving an early max after an adjacent later sum can only decrease the nested scalar. Repeated swaps reach a sum-before-max order. Thus any mixed order gives an upper bound on the true marginal-MAP mass; the all-sum-first schedule is sufficient for exactness.

The flag `sum_before_max` is deliberately sufficient, not necessary. Disconnected factors, exact numerical ties and other structure can make interleavings exact. `matches_exact_value` is an observed comparison against the exhaustive reference, not a theorem for unseen factors. The code does not recover a feasible MAP assignment from an invalid schedule. The normalized elimination scalar is still a bound and can exceed one; it is not labeled a probability.

For `[6,5,0,5]` on `[0,1]`, x0 MAP and x1 hidden, the correct masses are `[6,10]`. Max-first produces 11. Joint MAP instead selects `(0,0)` with weight 6. This original control separates three questions without attributing a learned mechanism.

Zero Z_e means no posterior exists. Weight maximizers are then algebraic zero ties, not meaningful posterior predictions. Probabilities and normalized bound are null; no epsilon repair or fallback normalization is used.

## Fill graph and work census

The primal graph has an edge for any pair of active variables in the same factor scope. Removing v connects its remaining neighbors into a clique; previously absent edges are fill edges. Structural scopes are kept even when a factor is all zero or numerically constant. Induced width is the largest neighbor count along this order. Enumerating all orders gives the exact minimum induced width for this tiny evidence-sliced primal graph, not a scalable treewidth algorithm or the constrained marginal-MAP width.

Joined entries at a step are `2^|join_scope|`; output entries are `2^(|join_scope|-1)`. Total joined entries sum those allocations. Peak joined/output entries describe the largest single allocated table, not total simultaneously live storage. Conceptual arithmetic counts use a minimal fold: m input factors require m-1 multiplications per joined entry, each binary sum one addition per output and each binary max one comparison per output. A final scalar fold takes m-1 multiplications. The implementation may multiply by an identity or compare an initial value; these counters are a transparent arithmetic model, not measured executed instructions. Counts exclude slicing, exact-rational bit growth, copying, Python/DOM overhead, all-order census and dense reference enumeration.

A four-leaf star on five variables gives 62 joined entries at width 4 for hub-first versus 18 at width 1 for leaves-first, with identical all-sum partition. The browser can replay any census order. No zero sparsity, heuristic ordering, runtime improvement, memory benchmark or scalability claim is made.

## Verification

The 24 Python tests use a direct dense assignment dictionary independent of factor helpers. It is repeatedly reduced over full assignments to check mixed-order results and every partial product of current factors. A separate adjacency-set algorithm checks every census width and joined allocation count. All source-factor indices and joined products, output branches, all max ties, evidence and units are checked. Forty scalar-rescaling controls preserve normalized probabilities and scale masses/gaps by 3/2. Invalid CLI input leaves an existing output untouched.

24 retained reports include the ordering counterexample, all-sum/all-max, two star schedules and evidence, zero mass, isolated variables, scalar/disconnected factors, all-observed inputs, exact ties and ten seeded cyclic graphs. Independent Python Fraction and JavaScript BigInt reports match all 23,114 exact leaves. Actual Chrome checks all 1,939 selected records (325 joint + 406 joined + 203 reductions + 75 steps + 69 sliced factors + 861 census records), 21 fill edges, 75 bars, tables/summaries, replay and exports. Responsive containment is checked at 320,390,1280; actual screenshots are inspected independently. Wheel tests reproduce every fixture and CLI from an isolated environment outside checkout and verify package origins.

## Established theory and sources

All sources below were read for definitions or scope, not copied as implementation. Research checked six primary/survey/course sources. Dates are those visible in the documents; no claim that this method is new.

1. Peyrard et al., *Exact or approximate inference in graphical models: why the choice is dictated by the treewidth, and how variable elimination can be exploited*, arXiv:1506.08544v2. Sections 2.2,3.2-3.3 support combination/elimination, evidence restriction, fill edges and induced width. [PDF](https://arxiv.org/pdf/1506.08544v2.pdf). The displayed prose example repeats a potential subscript; our generic union-of-scopes implementation does not import that typography. The paper's efficient-inference discussion is not a runtime claim about this exhaustive inspector.
2. Qiang Liu and Alexander Ihler, *Variational Algorithms for Marginal MAP*, arXiv:1202.3742. Introduction and section 3 define mixed inference and noncommuting sum/max schedules. [PDF](https://arxiv.org/pdf/1202.3742). Their variational/TRW/CCCP algorithms and novel bounds are not implemented or claimed here. Our elementary sum/max inequality is documented above, not attributed as their variational result.
3. Raquel Urtasun and Tamir Hazan, *Probabilistic Graphical Models*, April 20,2011 lecture. Evidence restriction and factor-allocation costs. [PDF](https://www.cs.toronto.edu/~urtasun/courses/GraphicalModels/lecture9.pdf). Its P(e) normalizer is for normalized probabilistic models; our arbitrary potentials require the Z_e qualification above. Approximate asymptotic slide counts are not copied into exact counters.
4. MIT 6.438, *Recitation 4: Elimination algorithm, reconstituted graph, triangulation*, Fall 2014. Given-potential queries, neighbor cliques and reconstituted graphs. [PDF](https://ocw.mit.edu/courses/6-438-algorithms-for-inference-fall-2014/ac39d786c86c4909bd83ce0de592375b_MIT6_438F14_rec4.pdf). No figure or pseudocode is copied.
5. Matt Gormley, *Factor Graphs and Exact Inference*, lecture22, November13,2016. Factor replacement, sum-product/max-product semirings and the distinction between max values and backtracked assignments. [PDF](https://www.cs.cmu.edu/~mgormley/courses/10601b-f16/lectureSlides/lecture22-pgm2.pdf). No tree-belief-propagation procedure is implemented.
6. Junkyu Lee, Radu Marinescu, Rina Dechter and Alexander Ihler, *From Exact to Anytime Solutions for Marginal MAP*,2016. Introduction/background define nonnegative-potential mixed queries and constrained orders. [PDF](https://ics.uci.edu/~dechter/publications/r225.pdf). The source states positive potentials; our zero-potential extension uses the same distributivity with explicit undefined-probability handling. No anytime/AND-OR/mini-bucket search or benchmark result is imported.

Implementation, independent tests and original synthetic controls by Midya Rahmani. No external source code, plots, data or weights copied. Not formal verification, probabilistic model learning, calibrated prediction or a new theorem.
