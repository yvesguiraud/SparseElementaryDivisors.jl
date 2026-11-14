# Documentation site of SparseElementaryDivisors.jl, built with Documenter.jl.
#
# Local build, from the root of the repository:
#     julia --project=docs -e "using Pkg; Pkg.instantiate()"
#     julia --project=docs docs/make.jl
# then open docs/build/index.html. Outside the CI, deploydocs publishes nothing; the site is
# published to GitHub Pages by the workflow .github/workflows/Documentation.yml.

using Documenter
using SparseElementaryDivisors

makedocs(;
    sitename = "SparseElementaryDivisors.jl",
    authors = "Yves Guiraud",
    modules = [SparseElementaryDivisors],
    # Every docstring of an exported or public name must be in the manual; the docstrings of
    # the internal functions stay in the code. A missing docstring or a broken cross-reference
    # makes the build fail (warnonly = false, the default).
    checkdocs = :public,
    doctest = false,
    format = Documenter.HTML(;
        # Pretty URLs (user_guide/ instead of user_guide.html) on the CI only, so that the
        # local build can be read by opening docs/build/index.html, without a web server.
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://yvesguiraud.github.io/SparseElementaryDivisors.jl/stable/",
        edit_link = "main",
        # Hides the headings of the current page in the sidebar (docs/src/assets/custom.css).
        assets = ["assets/custom.css"],
    ),
    pages = [
        "Home" => "index.md",
        "User guide" => "user_guide.md",
        "Benchmarks" => "benchmarks.md",
        "Reference manual" => "reference_manual.md",
        "Bibliography" => "bibliography.md",
    ],
)

# main gives dev/, a tag vX.Y.Z gives vX.Y.Z/ and stable/; nothing is published for a pull
# request. Authentication: the deploy key DOCUMENTER_KEY if the secret is set, otherwise
# GITHUB_TOKEN.
deploydocs(;
    repo = "github.com/yvesguiraud/SparseElementaryDivisors.jl.git",
    devbranch = "main",
)
