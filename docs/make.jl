using Documenter, DocumenterVitepress

using Literate

# build Quarto notebooks

for folder in ("Raster-in-45-min",)
    for qmd_file in filter(endswith(".qmd"), readdir(joinpath(dirname(@__DIR__), folder); join = true))
        dir, file = splitdir(qmd_file)
        cd(dir) do
            run(`quarto render $(qmd_file)`)
        end
        md_file = replace(qmd_file, ".qmd" => ".md")
        cp(
            md_file, 
            joinpath(@__DIR__, "src", "content", folder, splitdir(md_file)[2]); 
            force = true
        )
    end
end