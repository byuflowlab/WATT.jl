using Documenter, Rotors

makedocs(sitename="Rotors.jl", pages = [
    "Getting Started" => "gettingstarted.md",
    "Examples" => ["riso_example.md"],
    "Developers" => "developers.md",
    "API Reference" => "apireference.md"
])

deploydocs(
    repo = "github.com/byuflowlab/Rotors.jl.git",
)