# FactorLab

Exact finite factor elimination with every product, reduction and fill edge. Dependency-free Python plus an independent BigInt JavaScript inspector.

[Live inspector](https://midya8.github.io/factorlab/) · [Method and sources](METHOD.md)

## Sum is not max

A two-variable factor on scope `[0,1]` has little-endian values `[6,5,0,5]`. Let x0 be the MAP decision and x1 be hidden:

- Sum x1 first: decision masses are `[6,10]`, so x0=1 is optimal with probability `10/16 = 5/8`.
- Maximize x0 first: `max(6,5) + max(0,5) = 11`. This is an upper bound, not an attainable marginal-MAP mass.
- The single most likely joint assignment is x0=0,x1=0 with weight 6. Projecting that joint mode picks the wrong marginal decision.

A five-variable star also separates correctness from work: hub-first and leaves-first all-sum elimination agree exactly, but allocate 62 versus 18 joined table entries and have induced width 4 versus 1. These are dense entry counts, not timing or memory benchmarks.

## Run offline

```sh
python -m unittest discover -s tests -v
python -m factorlab examples/demo.json --output witness.json
python tools/build_audits.py
npm test
bash tools/wheel_smoke.sh
python -m http.server --directory docs 8000
```

Python 3.10+ has no runtime dependencies. The Node tests use Node 22+ and an installed `google-chrome`, or `CHROME_BIN`. The normal inspector runs without a model/API call; browser controls are fetched from the same static site.

## What stays inspectable

- Evidence slicing, scalar factors and explicit unit factors for absent variables.
- Every joined assignment and source-factor local index/value/product.
- Every two-branch sum or max reduction and all maximizing ties.
- Each elimination scope, new fill edge, remaining graph and dense allocation count.
- Independent complete enumeration, partition, joint-mode ties and marginal-MAP assignment masses.
- Every small elimination order, with a sufficient sum-before-max flag, induced width and entry counts.
- Editable inputs, order replay, current-only JSON export and stale/error clearing.

## Verification

24 Python tests on 3.10-3.13 include 160 direct dense scalar oracles, 100 full partial-factor-product oracles, every order on 20 random four-variable graphs for all-sum/all-max/mixed queries, 40 positive-scale controls, complete join/reduction checks and independent graph-degree census checks. Controls within a test are not separate test methods.

24 complete native/BigInt reports match 23,114 exact leaves. Actual Chrome checks 325 joint records, 406 joined products, 203 reductions, 75 step records, 69 sliced factors, 861 census records, 21 fill edges and 75 bars. It also checks all summary/table values, SVG edge coordinates/colors, order replay, export, stale/invalid behavior and 320/390/1280 containment. Installed-wheel reproduction checks every report, CLI and module origin outside checkout.

## Boundaries

Five binary variables, ten nonnegative rational factors, up to 32 complete assignments and 120 orders. Exhaustive references and order census are exponential. This is established inference on original synthetic potentials, not model learning, prediction quality, scalable inference, a new algorithm, belief propagation, variational inference or anytime search. The mixed-order scalar is an upper bound, not an assignment. Zero mass leaves probabilities null. Details, exact input conventions, conceptual count definitions and prior-work qualifications are in [METHOD.md](METHOD.md).

Implementation, independent checks and original controls by Midya Rahmani. No external code, figures, datasets or weights copied.
