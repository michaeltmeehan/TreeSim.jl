# Validation Plan

## Purpose

This document defines the validation work required for `TreeSim.jl` to be considered a trustworthy standalone tree package within the recovered outbreak-modelling ecosystem.

`TreeSim.jl` is the canonical package for tree representation, traversal, and validation. Its responsibility is to provide a minimal, stable, and semantically clear tree core that downstream packages may rely upon.

---

## Current scope

`TreeSim.jl` currently supports:

- a canonical `Tree` type
- node typing via `NodeKind`
- core traversal utilities
- tree validation utilities
- minimal stable standalone operation independent of the broader ecosystem

---

## Out of scope

The following are explicitly out of scope for the current phase:

- birth-death likelihood calculations
- epidemic simulation
- sequence simulation
- plotting and visualization layers
- heavy conversion/import/export frameworks
- orchestration logic

---

## Core validation questions

Validation work must establish that:

1. a tree that passes validation is structurally coherent
2. tree semantics are sufficiently explicit for downstream consumers
3. traversal utilities behave correctly across all supported topologies
4. malformed or inconsistent trees are rejected reliably
5. `NodeKind` semantics are stable and unambiguous

---

## Required validation areas

### 1. Structural invariants

Validation tests must cover:

- exactly one root where required
- valid parent-child relationships
- no impossible parent references
- no disconnected components in valid trees
- no duplicate or conflicting ancestry
- consistency of left/right child fields
- consistency of parent pointers

### 2. Temporal invariants

Time-ordered nodes are a required invariant of the canonical tree representation.

Validation tests must therefore cover:

- valid temporal ordering between every parent and child pair
- rejection of any tree in which a child occurs earlier than its parent under the package's time convention
- explicit handling of whether tied parent-child times are permitted or forbidden
- rejection of any tree that violates the package's canonical node ordering convention if nodes are also required to appear in temporal order in storage

This invariant is intentionally strict. It is imposed not only for structural clarity, but also to simplify downstream traversal, validation, and analytical processing.

### 3. Node kind semantics

Tests must verify:

- permitted structural roles of each `NodeKind`
- forbidden structural roles of each `NodeKind`
- compatibility between node type and child configuration
- consistency of sampled/internal/unary/binary interpretations

### 4. Traversal correctness

Tests must verify:

- preorder/postorder/reverse traversal correctness where supported
- root-based traversal behaviour
- coverage of all nodes exactly once where appropriate
- stable behaviour on degenerate but valid trees

### 5. Validation failure behaviour

Tests must verify that malformed trees fail clearly for cases such as:

- missing root
- multiple roots where unsupported
- cyclic ancestry
- invalid child indices
- invalid parent indices
- impossible node kinds
- temporal inconsistency
- inconsistent child-parent cross-references

---

## Canonical fixture suite

A validation fixture suite should include hand-constructed examples for:

- minimal valid tree
- simple binary tree
- unary internal tree
- sampled leaf example
- malformed parent reference
- malformed child reference
- disconnected tree
- time-inconsistent tree
- node-kind inconsistent tree
- impossible root configuration

These fixtures should be reused throughout the test suite.

---

## Downstream validation relevance

Because `BDUtils.jl` depends on `TreeSim.jl`, validation must also establish a clear contract for what constitutes an admissible tree for downstream analytical use.

This package does **not** validate birth-death admissibility itself unless explicitly implemented, but it must provide the structural and semantic guarantees that downstream consumers rely on.

---

## Exit criteria for Phase 2

`TreeSim.jl` should not be considered fully validated for this phase until:

- tree semantics are documented clearly
- all supported `NodeKind` values have explicit structural meaning
- a canonical fixture suite exists
- traversal behaviour is thoroughly tested
- malformed trees are rejected reliably
- downstream object assumptions are written down
- the public API for the minimal tree core is stable

---

## Evidence of successful validation

Evidence should include:

- unit tests for canonical valid trees
- unit tests for canonical invalid trees
- traversal tests with exact expected outputs
- validation tests covering each invariant category
- documentation aligning code behaviour with the stated semantics
