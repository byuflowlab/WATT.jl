Create a focused implementation sub-plan for a single phase of the WATT.jl rework.

## Instructions

The user will specify a phase number (1–8). Your job is to produce a detailed, session-scoped sub-plan that can be executed start-to-finish in one working session.

### Step 1: Read the long-term plan
Read `plan.md` in the project root. Identify the target phase and its goals, entry/exit criteria, tasks, and dependencies.

### Step 2: Explore relevant source files
Read the source files listed in the phase's "Critical Files" or task list. Focus on:
- The exact current state of the code (what exists, what is dead, what is broken)
- Function signatures, struct fields, and call sites that will be affected
- Any existing tests that cover the area
- Any patterns already established in the codebase that should be followed

Also read:
- `.claude/context/GXBEAM_STYLE.md` — for style/API conventions to apply
- `.claude/context/JULIA_GUIDE.md` — for Julia performance and style rules

### Step 3: Write the sub-plan

The sub-plan must include:

**Header**
- Phase number and title
- Status: In Progress
- Date
- Link back to `plan.md`

**Entry check**
- Confirm the phase's entry criteria are met (or note what is blocking)
- List any unresolved questions to ask the user before starting

**Ordered task list**
- Number each task
- Mark which tasks are independent and can be parallelized vs. must be sequential
- For each task: file path(s), what to change, what to watch out for (AD compatibility, type stability, naming)
- Flag any task that touches a known active bug from the plan

**Test plan**
- What tests to run after each task group to verify correctness
- For Phase 1: run the golden-value script (user must provide) and verify outputs match reference
- For Phases 2+: list specific test files and what assertions to check
- Note if AD tests should be run (only when `ENV["WATT_AD_TESTS"] == "true"`)

**Exit checklist**
- Checklist derived from the phase's exit criteria in `plan.md`
- Include: tests pass, no regressions vs. golden reference, style matches GXBEAM_STYLE.md

**Next session note**
- One sentence on what Phase comes next and what the entry criteria will be

### Output format

Write the sub-plan as a markdown document. Do not execute any of the tasks — only produce the plan. The user will approve it before implementation begins.
