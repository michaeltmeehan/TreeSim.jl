# Trust Criteria

## Role of this package

`TreeSim.jl` is the canonical standalone tree package in the recovered outbreak-modelling ecosystem.

Its role is to provide a minimal and trustworthy tree representation layer for downstream analytical and simulation packages.

---

## Trust goal

`TreeSim.jl` is trustworthy when a downstream package can safely assume that a tree passing validation is:

- structurally coherent
- semantically interpretable
- traversable without ambiguity
- stable with respect to the documented minimal public API

---

## What trust does mean here

Trust in this phase means:

- the package loads cleanly
- the minimal public tree API is intentional
- tree structure invariants are explicit and tested
- `NodeKind` semantics are documented and enforced
- traversal methods behave correctly on supported trees
- malformed trees are rejected rather than silently tolerated

---

## What trust does not mean here

Trust in this phase does **not** imply:

- support for all conceivable tree encodings
- built-in visualization or plotting
- built-in likelihood calculations
- ecosystem-wide integration
- support for rich import/export layers
- full phylogenetic or epidemiological semantics beyond the current minimal tree core

---

## Stable trust boundary

The current trust boundary includes:

- `Tree`
- `NodeKind`
- traversal utilities included in the active core
- validation utilities included in the active core

Anything outside this boundary should be considered provisional unless explicitly promoted into the stable API.

---

## Conditions required for trust

### 1. Semantic clarity

The meaning of all node kinds and structural fields must be explicitly documented.

### 2. Structural reliability

Valid trees must satisfy the documented invariants, and invalid trees must fail validation clearly.

### 3. Traversal reliability

Traversal utilities must return correct and reproducible outputs across all supported tree shapes.

### 4. API stability

The minimal public interface must be intentional and protected from casual expansion.

### 5. Downstream safety

A downstream consumer should be able to rely on validated trees without needing to reverse-engineer internal assumptions.

---

## Known trust risks

Current or likely risks include:

- under-specified `NodeKind` semantics
- edge cases involving unary or sampled-node structures
- traversal mistakes on unusual topologies
- malformed trees that are superficially plausible
- ambiguity between structural and downstream analytical validity

---

## Required evidence before calling this package trustworthy

The following evidence is required:

- a written semantics note for core tree objects
- canonical valid and invalid fixture trees
- comprehensive invariant tests
- traversal tests with expected outputs
- explicit downstream contract notes where relevant

---

## Phase-2 completion standard

For the purposes of this project phase, `TreeSim.jl` is trustworthy when:

1. its minimal stable core is clearly scoped
2. its semantics are written down
3. its core invariants are enforced by tests
4. malformed structures are reliably rejected
5. downstream packages can safely treat it as the canonical tree layer
