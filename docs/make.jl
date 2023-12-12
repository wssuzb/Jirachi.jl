using Jirachi
using Documenter

DocMeta.setdocmeta!(Jirachi, :DocTestSetup, :(using Jirachi); recursive=true)

makedocs(;
    modules=[Jirachi],
    authors="zhenbo su",
    repo="https://github.com/wssuzb/Jirachi.jl/blob/{commit}{path}#{line}",
    sitename="Jirachi.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wssuzb.github.io/Jirachi.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/wssuzb/Jirachi.jl",
    devbranch="main",
)
