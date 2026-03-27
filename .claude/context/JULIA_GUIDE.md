# Julia Style, Performance & Scoping Reference

> Condensed from the official Julia documentation:
> - https://docs.julialang.org/en/v1/manual/performance-tips/
> - https://docs.julialang.org/en/v1/manual/style-guide/
> - https://docs.julialang.org/en/v1/manual/variables-and-scoping/

---

## Style

### Naming Conventions
- **Modules and types**: `UpperCamelCase` — `module SparseArrays`, `struct UnitRange`
- **Functions**: `lowercase`, words squashed together where readable — `isequal`, `haskey`; underscores used for combinations — `remotecall_fetch`
- **Mutating functions**: append `!` — `sort!`, `push!`, `pop!`. Exception: IO and RNG functions mutate by nature, so `!` signals mutation *beyond* advancing IO/RNG state (e.g. `rand!(x)` mutates both the RNG and `x`)
- **"Internal" names**: prefix or suffix with `_`, though this is convention not a rule

### Function Design
- **Write functions, not scripts.** Top-level scripted code is slow and hard to reuse. Wrap everything in functions.
- **Functions take arguments** rather than operating on globals (except true constants like `pi`).
- **Break compound functions into multiple dispatch definitions** rather than using `isa` branches:
  ```julia
  # Bad
  function mynorm(A)
      if isa(A, Vector); ...; elseif isa(A, Matrix); ...; end
  end

  # Good
  mynorm(x::Vector) = sqrt(real(dot(x, x)))
  mynorm(A::Matrix) = maximum(svdvals(A))
  ```
- **Be generic with argument types.** Prefer `Integer` over `Int`, `AbstractArray` over `Array`. Often no annotation at all is best — Julia dispatches on concrete types automatically with no performance penalty.
  ```julia
  addone(x) = x + oneunit(x)   # works for any type supporting + and oneunit
  ```
- **Handle type diversity in the caller**, not inside the function. If `foo` truly needs `Int`, declare `foo(x::Int, y::Int)` and let the caller do `foo(Int(x), Int(y))`.
- **Argument ordering** (follow Julia Base convention):
  1. Function argument (enables `do` blocks)
  2. IO stream
  3. Mutated input
  4. Type
  5. Non-mutated input
  6. Key / index
  7. Value
  8. Varargs
  9. Keyword arguments

### Type Usage
- **Avoid overly-specific types** in function signatures. Use abstract types or no annotation.
- **Avoid strange `Union` types** like `Union{Function,AbstractString}` — usually a sign of a design problem.
- **Avoid elaborate container types** like `Vector{Union{Int,AbstractString,Tuple,Array}}`. Prefer `Vector{Any}` and annotate at use sites instead.
- **Prefer exported methods over direct field access.** Treat a module's exported functions as its public API. Direct field access couples you to implementation details.
- **Use `isa` and `<:` for type testing**, not `==`. Exact type equality (`T == Float64`) is only appropriate for known concrete types.
- **Avoid unnecessary static parameters:**
  ```julia
  # Unnecessary
  foo(x::T) where {T<:Real} = ...

  # Better (identical performance)
  foo(x::Real) = ...
  ```
- **Constructors should return an instance of their own type.** A `T(x)` call is expected to return a `T`.

### Other Conventions
- **Indentation**: 4 spaces per level.
- **No parens around `if`/`while` conditions**: `if a == b`, not `if (a == b)`.
- **Don't wrap named functions in trivial lambdas**: `map(f, a)`, not `map(x->f(x), a)`.
- **Don't overuse `...` (splat)**: `[a; b]` concatenates arrays directly; prefer `collect(a)` over `[a...]`.
- **Don't overuse macros** — ask whether a function would do. Avoid `eval` inside macros.
- **Avoid type piracy** — don't extend methods in Base or other packages on types you didn't define.
- **Don't overload methods of base container types** like `Vector` — users expect predictable behavior.
- **Avoid floats in generic numeric code**:
  ```julia
  f(x) = 2.0 * x   # promotes Int → Float64 unexpectedly
  g(x) = 2 * x     # preserves type via promotion rules
  ```
- **Don't overuse `try`/`catch`** — avoid errors rather than catching them.
- **Mark unsafe operations** with `unsafe` in the name or add runtime safety checks before exposing them.

---

## Performance

### The Golden Rules

**1. Performance-critical code must be inside a function.**
Julia's compiler heavily optimizes function-scoped code. Top-level code is treated as dynamic global scope and cannot be optimized the same way.

**2. Avoid untyped global variables.**
Untyped globals can change type at any time, preventing compiler optimization.
```julia
# Bad: compiler can't infer type of x
x = rand(1000)
function sum_global()
    s = 0.0
    for i in x; s += i; end
end

# Good: pass as argument
function sum_arg(x)
    s = 0.0
    for i in x; s += i; end
end

# Also good: declare constant or typed global
const WEIGHTS = rand(1000)
global typed_x::Vector{Float64} = rand(1000)
```

**3. Unexpected memory allocation signals a problem.**
Use `@time` (after a warm-up call) and watch for allocations. Zero allocations in a numeric loop is achievable and expected. Use `BenchmarkTools.jl` for rigorous measurement.

### Type Stability

**Write type-stable functions** — functions where the return type can be inferred from argument types alone at compile time.

**Avoid containers with abstract type parameters:**
```julia
Real[]         # bad: array of pointers, slow
Float64[]      # good: contiguous memory, fast
```

**Avoid abstract-typed struct fields:**
```julia
# Bad: compiler can't specialize on `a`
struct MyType
    a         # inferred as Any
end

# Good: concrete type
struct MyType
    a::Float64
end

# Also good: parametric for flexibility without sacrificing performance
struct MyType{T}
    a::T
end
```

**Avoid changing a variable's type** within a function body — this defeats type inference.

**Use function barriers** to isolate type-unstable code. If you must work with a type-unstable container, extract the inner loop into a separate function so the compiler can specialize on the concrete element type it receives.

**Annotate values from untyped locations** at the point of use:
```julia
for i in x::Vector{Float64}   # tells compiler what x contains
```

**Use `@code_warntype`** to spot type instability — red/yellow annotations indicate inference failures.

### Memory & Arrays

**Pre-allocate outputs** rather than growing arrays in a loop:
```julia
# Bad: allocates on every iteration
result = []
for x in data
    push!(result, f(x))
end

# Good: pre-allocate
result = Vector{Float64}(undef, length(data))
for (i, x) in enumerate(data)
    result[i] = f(x)
end
# Or simply:
result = f.(data)
```

**Use views for array slices** to avoid copying:
```julia
sum(view(A, :, 3))   # no copy
# or
@views sum(A[:, 3])
```

**Access arrays in column-major order** (Julia is column-major like Fortran, unlike C/Python):
```julia
# Bad: row-major traversal
for i in 1:n, j in 1:m
    A[i,j] = ...
end

# Good: column-major traversal
for j in 1:m, i in 1:n
    A[i,j] = ...
end
```

**Fuse broadcasts with dot syntax** — Julia fuses adjacent broadcasts into a single loop:
```julia
y .= a .* x .+ b   # one loop, no temporaries
```

**Consider `StaticArrays.jl`** for small, fixed-size arrays (e.g. 3D vectors, 4×4 matrices) — avoids heap allocation entirely.

### Performance Annotations
- `@inbounds` — disable bounds checking in tight loops (only after verifying correctness)
- `@fastmath` — allow floating-point reassociation (may change results slightly)
- `@simd` — hint for SIMD vectorization on innermost loops
- `@views` — convert all slices in an expression to views

### Profiling Tools
| Tool | Purpose |
|------|---------|
| `@time` | Quick timing + allocation count (run twice; first call includes compilation) |
| `BenchmarkTools.@benchmark` | Rigorous statistical benchmarking |
| `Profile.@profile` + `ProfileView.jl` | Sampling profiler with flame graph |
| `@code_warntype` | Spot type instability (look for red/yellow) |
| `@code_native` / `@code_llvm` | Inspect generated machine code |
| `JET.jl` | Static analysis for type errors and performance issues |
| `--track-allocation=user` | Find exact lines causing heap allocation |

---

## Scoping

### Scope Block Summary

| Construct | Scope type | Notes |
|-----------|-----------|-------|
| `module`, `baremodule` | Global | Each module is its own namespace |
| `function`, `do`, `let`, comprehensions, generators | Local (hard) | Assignment always creates a new local |
| `struct`, `macro` | Local (hard) | |
| `for`, `while`, `try`/`catch` | Local (soft) | Assignment may interact with enclosing globals |
| `begin`, `if` | *None* | Do NOT introduce a new scope |

Julia uses **lexical scoping** — a function's scope is determined by where it is *defined*, not where it is *called*.

### The Core Rules

**Hard scope** (function body, `let`, comprehension, generator, `struct`, `macro`): assignment to `x` always creates a new local `x`, regardless of any global `x`.

**Soft scope** (`for`, `while`, `try`/`catch` outside a function): behavior depends on context:
- In the **REPL/interactive**: if a global `x` exists, assignment updates it; otherwise creates a new local.
- In **files/scripts**: if a global `x` exists, Julia emits an ambiguity warning and creates a new local anyway — leading to `UndefVarError` when the local is used before assignment.

**The practical consequence** — this script-level pattern is a common bug:
```julia
# In a file — BROKEN
s = 0
for i = 1:10
    s += i    # Warning: ambiguous — creates new local `s`, not global
end
println(s)    # prints 0, not 55 — global `s` was never updated
```

Fix it by wrapping in a function (preferred) or using `global` explicitly:
```julia
# Fix 1: wrap in a function (preferred)
function sum_to(n)
    s = 0
    for i = 1:n
        s += i
    end
    return s
end

# Fix 2: explicit global (acceptable for scripts)
s = 0
for i = 1:10
    global s += i
end
```

### Key Behaviors to Remember

- `for` loop variables (`i` in `for i = ...`) are **local to the loop body** — they don't exist after the loop.
- Variables assigned inside a `for` loop body that are **not** already defined in an outer local scope create new loop-local variables.
- `begin`/`if` blocks **do not create a new scope** — variables assigned inside are visible outside.
- Inner functions **close over** outer local variables and can mutate them:
  ```julia
  function outer()
      x = 0
      inner() = (x += 1)   # closes over x
      inner(); inner()
      return x              # returns 2
  end
  ```
- The entire enclosing local scope is parsed before inner local meanings are resolved — a variable's scope is determined statically, not in execution order.

### Constants and Typed Globals
```julia
# Constant — type is inferred, value never changes, fully optimized
const MAX_ITER = 1000
const LOOKUP = Dict("a" => 1, "b" => 2)

# Typed global — value can change, but type is fixed, enables optimization
global counter::Int = 0
```

`const` globals do not need type annotations (type is inferred). Non-constant globals should be typed if used in performance-sensitive code, or better yet, passed as function arguments.

---

## Quick Checklist

Before committing Julia code, verify:

- [ ] Performance-critical code is inside a function
- [ ] No untyped globals used in hot paths (use `const`, typed globals, or function arguments)
- [ ] Struct fields have concrete or parametric types, not abstract types
- [ ] Arrays are typed concretely (`Float64[]`, not `Real[]`)
- [ ] Functions are type-stable (`@code_warntype` shows no red)
- [ ] Array traversal is column-major (outer loop = columns, inner loop = rows)
- [ ] Slices use `@views` where allocation matters
- [ ] Broadcasts are fused with dot syntax
- [ ] Mutating functions end in `!`
- [ ] No ambiguous soft-scope assignments in script-level code
- [ ] No type piracy
