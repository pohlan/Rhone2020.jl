#### HOW TO CREATE MARKDOWN FILES ####

using Literate
cd(@__DIR__)

# execute=true produces plots

#Literate.markdown("paper_figures.jl"; execute=true, documenter=false)
#Literate.markdown("other_figures.jl"; execute=true, documenter=false)
#Literate.markdown("raw_data.jl"; execute=true, documenter=false)
#Literate.markdown("calibration_plots.jl"; execute=true, documenter=false)
#Literate.markdown("derived_quantities.jl"; execute=true, documenter=false)


# to convert .md files to pdfs, linux command, e.g.:
# $ pandoc -s -V geometry:margin=1in -o derived_quantities.pdf derived_quantities.md
