using Documenter, PeriodicMatrices, LinearAlgebra
DocMeta.setdocmeta!(PeriodicMatrices, :DocTestSetup, :(using PeriodicMatrices); recursive=true)

makedocs(warnonly = true, 
  modules  = [PeriodicMatrices],
  sitename = "PeriodicMatrices.jl",
  authors  = "Andreas Varga",
  format   = Documenter.HTML(prettyurls = false),
  pages    = [
      "Home"   => "index.md",
      "Library" => [ 
         "Data Types and Constructors" => [
          "pmtypes.md"
          ],
          "Periodic Matrix Operations" => "pmops.md",
          "Periodic Matrix Tools" => "pmtools.md",
          "Periodic matrix conversions" => "pmconv.md"
         ],
         "Internals" => [
         "Utilities" => "pmutilities.md",
         "Periodic Schur Decomposition" => "pschur.md",
         "SLICOT Based Wrappers" => "slicot.md"
         ],
     "Index" => "makeindex.md"
  ]
)

deploydocs(
  repo = "github.com/andreasvarga/PeriodicMatrices.jl.git",
  target = "build",
  devbranch = "main"
)
