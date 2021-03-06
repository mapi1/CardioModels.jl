using Documenter, CardioModels

makedocs(
    modules = [CardioModels],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Marius Pille",
    sitename = "CardioModels.jl",
    pages = [
        "General" => "index.md",
        "DeBoer" => "DeBoer.md",
        "Karemaker" => "Karemaker.md",
        "Baselli" => "Baselli.md",
    ]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/mapi1/CardioModels.jl.git",
    push_preview = true
)
