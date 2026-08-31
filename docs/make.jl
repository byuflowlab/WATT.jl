using Documenter, WATT

makedocs(
    sitename = "WATT.jl",
    modules = [WATT],
    # Fail the build on a missing docstring rather than shipping a hole in the
    # API reference. Every exported symbol must appear in an @docs block.
    checkdocs = :exports,
    pages = [
        "Getting Started" => "gettingstarted.md",
        "Examples" => [
            "aeroonly.md",
            "steady.md",
            "gradients.md",
            "assembly.md",
        ],
        "Developers" => "developers.md",
        "API Reference" => ["apireference.md", "gpu.md"],
    ],
)

deploydocs(
    repo = "github.com/byuflowlab/WATT.jl.git",
    devbranch = "master",
    # Publish a preview under previews/PR## for each pull request. The
    # cleanup.yaml workflow already expects these to exist and deletes them
    # when the PR closes.
    push_preview = true,
)
