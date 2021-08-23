#### HOW TO CREATE MARKDOWN FILES ####

using Literate, PyPlot
cd(@__DIR__)
pygui(false) # no interactive plots

tmp_outputdir = "tmp"
mkpath(tmp_outputdir)
pdfdir = joinpath(@__DIR__, "../../../products/additional_figures/")
mkpath(pdfdir)

# Run Literate.jl over the files:
# (execute=true produces plots)
Literate.markdown("raw_data.jl", tmp_outputdir; execute=true, documenter=false)
Literate.markdown("calibration_plots.jl", tmp_outputdir; execute=true, documenter=false)
Literate.markdown("derived_quantities_extra.jl", tmp_outputdir; execute=true, documenter=false)
Literate.markdown("other_figures.jl", tmp_outputdir; execute=true, documenter=false)

# to convert .md files to pdfs, linux command, e.g.:
cd(tmp_outputdir)
for fl in readdir()
    if splitext(fl)[2]==".md"
        fll = splitext(fl)[1]
        try
            run(`pandoc -s -V geometry:margin=1in -o $fll.pdf $fll.md`)
            cp("$fll.pdf", joinpath(pdfdir, "$fll.pdf"), force=true)
        catch e
            println(e)
            println("Could not convert markdown-file $fll.md to pdf")
        end
    end
end

cd(@__DIR__)
#
rm(tmp_outputdir, recursive=true)
