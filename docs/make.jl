using Documenter, WATT

makedocs(sitename="WATT.jl", pages = [
    "Getting Started" => "gettingstarted.md",
    # "Examples" => ["steady.md"], #Todo: Is there a different way to tackle this. 
    "Developers" => "developers.md",
    "API Reference" => "apireference.md"
])

deploydocs(
    repo = "github.com/byuflowlab/WATT.jl.git",
    devbranch = "master",
)