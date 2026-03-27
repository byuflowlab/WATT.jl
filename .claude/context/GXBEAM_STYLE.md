# GXBeam.jl Style Reference

Patterns extracted from GXBeam.jl v0.3.1 source (`~/.julia/packages/GXBeam/525aD/src/`).
Use these as the target style for WATT.jl refactoring.

---

## Struct Conventions

### Parametric structs — always use `{TF}` for floating point type

```julia
struct PointState{TF}
    u::SVector{3, TF}
    theta::SVector{3, TF}
    F::SVector{3, TF}
    M::SVector{3, TF}
end

struct Element{TF}
    L::TF
    x::SVector{3, TF}
    compliance::SMatrix{6,6,TF,36}
    mass::SMatrix{6,6,TF,36}
end
```

- Type parameter is always named `TF` (floating point type)
- Additional container type parameters use `TV`, `TM`, `TP`, `TC`, `TE` etc.
- Fully-constrained signatures for complex structs:
  ```julia
  struct Assembly{TF, TP<:AbstractVector{<:AbstractVector{TF}},
      TC<:AbstractVector{<:Integer}, TE<:AbstractVector{Element{TF}}}
  ```
- Small fixed-size arrays: `SVector` / `SMatrix` from StaticArrays.jl (avoids heap allocation)
- Large/variable arrays: `AbstractVector{TF}` (allows both `Vector` and sparse/offset arrays)

### Immutable vs mutable

- **Immutable** (`struct`): pure data containers — state outputs, element properties, boundary conditions, load structs
- **Mutable** (`mutable struct`): system state that is updated in-place during solve — `StaticSystem`, `DynamicSystem`

The key distinction: if you only mutate array *contents* (not array references), immutable is fine. Use `mutable struct` only when field references themselves are reassigned.

### Always define `Base.eltype`

```julia
Base.eltype(::Element{TF}) where TF = TF
Base.eltype(::Type{Element{TF}}) where TF = TF   # both instance and type forms
```

### Type conversion constructor + `Base.convert`

```julia
Element{TF}(e::Element) where {TF} = Element{TF}(e.L, e.x, e.compliance, e.mass, e.Cab, e.mu)
Base.convert(::Type{Element{TF}}, e::Element) where {TF} = Element{TF}(e)
```

This enables `convert(Element{Float32}, elem)` and type-parameterized construction from existing instances — important for AD where dual numbers replace Float64.

### Outer constructors infer `TF` via `promote_type`

```julia
function Element(L, x, compliance, mass, Cab, mu)
    TF = promote_type(typeof(L), eltype(x), eltype(compliance), eltype(mass), eltype(Cab), eltype(mu))
    return Element{TF}(L, x, compliance, mass, Cab, mu)
end
```

Users call the outer constructor; the `{TF}` constructor is the canonical low-level form.

---

## Module Organization (`GXBeam.jl`)

### All exports live in the main module file — none in subfiles

```julia
module GXBeam

# ... using statements ...

export AbstractSystem, StaticSystem, DynamicSystem
export static_analysis, static_analysis!
export steady_state_analysis, steady_state_analysis!
# ... all exports grouped by category ...

# module-level constants
const GAUSS_NODES = SVector(-0.861..., ...)

# includes with one comment per group
include("math.jl")          # common math functions
include("assembly.jl")      # assembly creation
include("loads.jl")         # prescribed conditions, distributed loads, point masses
include("system.jl")        # system storage and pointers
include("analyses.jl")      # system analyses

# interface extensions in subdirectory
include("interfaces/diffeq.jl")
include("interfaces/reversediff.jl")
include("interfaces/writevtk.jl")

end
```

### Private/internal functions — no leading underscore, just don't export

GXBeam does not use `_foo` naming for private functions. Private functions simply aren't exported. No `@private` or similar macros.

### RecipesBase is a hard dependency in GXBeam

`using RecipesBase` appears directly in the main module. For WATT.jl we intentionally differ here: use `Requires.jl` conditional loading (see `plan.md` Phase 8).

---

## Dispatch Patterns

### Paired allocating / non-allocating functions

```julia
# Allocating (user-facing): creates system, calls mutating version
function static_analysis(assembly; kwargs...)
    system = StaticSystem(assembly)
    return static_analysis!(system, assembly; kwargs..., reset_state=false)
end

# Non-allocating (pre-allocated): takes pre-built system, mutates in place
function static_analysis!(system::StaticSystem, assembly; kwargs...)
    # ... solve ...
end
```

Pattern: `foo(args)` allocates and delegates to `foo!(preallocated, args)`. The `!` version is the workhorse; the non-`!` version is the convenience wrapper.

### Dispatch over abstract type, then concrete type

```julia
# First dispatch: accept any AbstractSystem
function update_body_acceleration_indices!(system::AbstractSystem, prescribed_conditions)
    update_body_acceleration_indices!(system.indices, prescribed_conditions)
    return system
end

# Second dispatch: work on the concrete SystemIndices
function update_body_acceleration_indices!(indices::SystemIndices, prescribed_conditions)
    # actual work here
    return indices
end
```

No `if isa(system, StaticSystem)` branching — dispatch handles it.

### Mutating functions return `self`

```julia
function update_body_acceleration_indices!(system::AbstractSystem, ...)
    # ... mutations ...
    return system   # ← always return the mutated object for chaining
end
```

---

## Docstring Format

GXBeam uses `# Section` headers (with `#`), not bold `**Section**`. For WATT.jl we use `**bold**` headers per user preference — note this divergence when writing docs.

### GXBeam actual format (for reference):

```julia
"""
    function_name(arg1, arg2; kwarg=default)

One-sentence or two-sentence description. May span lines.

# Arguments
 - `arg1`: Description (no type annotation in docstring)
 - `arg2`: Description

# Keyword Arguments
 - `kwarg = default`: Description

# Fields (for structs)
 - `field`: Description
"""
```

### WATT.jl adapted format (use this):

```julia
"""
    function_name(arg1, arg2; kwarg=default) -> out1, out2

One-sentence description.

**Arguments**
- `arg1::Type`: Description
- `arg2::Type`: Description

**Keyword Arguments**
- `kwarg::Type = default`: Description

**Returns**
- `out1::Type`: Description

**Notes**
AD compatibility notes, performance caveats, references.
"""
```

Key differences from GXBeam: bold headers, explicit return types in synopsis, type annotations on args.

### Struct docstring:

```julia
"""
    MyStruct{TF}

One-sentence description.

**Fields**
- `field1::SVector{3,TF}`: Description
- `field2::TF`: Description
"""
struct MyStruct{TF}
    ...
end
```

---

## Keyword Argument Patterns

### Default values use StaticArrays for zero vectors

```julia
function foo(assembly;
    gravity = (@SVector zeros(3)),          # ← not zeros(3)
    prescribed_conditions = Dict{Int, PrescribedConditions{Float64}}(),
    time = 0.0,
    reset_state = true,
    show_trace = false,
    method = :newton,
    ftol = 1e-9,
    iterations = 1000,
    )
```

### Optional parameters use `nothing` + `isnothing` guard

```julia
function Element(points, start, stop;
    frame = nothing,
    compliance = nothing,
    mass = nothing)

    if isnothing(compliance)
        compliance = @SMatrix zeros(6,6)
    end
    ...
end
```

---

## Other Notable Patterns

### `@unpack` for struct/NamedTuple destructuring

```julia
@unpack force_scaling, indices = system
```

Use `UnPack.jl` `@unpack` rather than manual field access when extracting multiple fields.

### NamedTuple packing with `(; ...)` syntax

```julia
constants = (;
    assembly, indices, two_dimensional, force_scaling,
    x=system.x, resid=system.r, converged=converged,
)
```

Semicolon-NamedTuple syntax for grouping closure arguments.

### `Ref` for scalar flags passed into closures/callbacks

```julia
converged = Ref(false)
# ... pass converged into callback ...
converged[] = true
```

### Abstract type hierarchy for system variants

```julia
abstract type AbstractSystem end
mutable struct StaticSystem{...}    <: AbstractSystem end
mutable struct DynamicSystem{...}   <: AbstractSystem end
mutable struct ExpandedSystem{...}  <: AbstractSystem end
```

Define an abstract supertype even for 2–3 concrete variants — enables dispatch without coupling call sites to concrete types.

### Interface extensions in `interfaces/` subdirectory

AD interfaces (ReverseDiff overloads, ChainRules), DifferentialEquations.jl integration, and WriteVTK integration live in `src/interfaces/`. This keeps the core physics files clean of third-party integration code.

---

## What GXBeam Does NOT Do

- No `isa` or `typeof` checks inside function bodies — all handled by dispatch
- No `@show` or `println` in any function
- No `using Revise` in tests
- No dead/commented code in source files
- Exported symbols are all in one place (`GXBeam.jl`), not scattered across files
- Does not use `NamedTuple` for primary data types — always a proper struct
