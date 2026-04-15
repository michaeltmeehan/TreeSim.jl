## `TreeSim.jl` next-phase plan

### Role

`TreeSim.jl` is the canonical standalone tree representation, traversal, and validation package.

### Current position

`TreeSim.jl` has been recovered and validated as a minimal trusted canonical tree core.

Tree simulation, especially broader historical coalescent-tree simulation recovery, is explicitly not a current priority.

### Immediate priorities

#### 1. Tree traversal utilities

- parent/child traversal helpers
- ancestor / descendant iterators
- preorder / postorder / breadth-first traversals
- leaf / internal-node iterators
- root-to-tip path utilities
- subtree extraction helpers

#### 2. Tree statistics

- number of leaves / internal nodes
- tree height
- branch-length summaries
- root-to-tip distance summaries
- balance-related summaries
- depth summaries
- subtree size summaries

#### 3. Visualization support

- simple tree plotting
- branch-length-aware plotting
- tip / node annotation support
- summary display helpers

#### 4. Processing utilities

- node filtering
- subtree extraction
- relabeling helpers
- annotation helpers
- tree simplification / pruning utilities

### Medium-term priorities

- tree mutation / perturbation utilities
- additional summary layers
- interchange / serialization helpers

### Explicit deferrals

- broad coalescent-tree simulation recovery
- turning `TreeSim.jl` into a general tree-simulation laboratory
- workflow orchestration logic

### Architectural guardrails

- preserve the identity of `TreeSim.jl` as the canonical tree package
- analytical enrichment is in scope
- broad generative expansion is not a current priority
