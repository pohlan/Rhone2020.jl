# Script which produces the overview map
#
# Requires command line tools:
# - unzip
# - gdal_merge.py
# used coordinate system: CH1903+/LV95

cd(@__DIR__)
using Rhone2020
# using Plots
# pyplot()
const R20 = Rhone2020

import PyPlot
const Plt = PyPlot
#Plt.pygui(true) # figure appears in window


using GeoData, Shapefile, ArchGDAL, DelimitedFiles, LaTeXStrings

# ------------------------------------------ #
#             Download data & unzip          #
#             not necessary anymore          #
# ------------------------------------------ #
#
# # surface DEM swisstopo SwissAlti3D (2m posting)
# files_surf = ["https://data.geo.admin.ch/ch.swisstopo.swissalti3d/swissalti3d_2019_267$(i)-116$(j)/swissalti3d_2019_267$(i)-116$(j)_2_2056_5728.tif" for i=1:3, j=0:2][:]
#
# # outline of Switzerland
# # https://www.swisstopo.admin.ch/de/geodata/landscape/boundaries3d.html
# zipfile_CH = "https://cms.geo.admin.ch/ogd/topography/swissBOUNDARIES3D.zip"
# folder_CH = "BOUNDARIES_2021/DATEN/swissBOUNDARIES3D/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.*"
#
# # Bed DEM and outline from SwissGlacierThickness-R2020
# # https://www.research-collection.ethz.ch/handle/20.500.11850/434697
# #
# # Rhone is SGI B43/03
# zipfile_outline = "https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/01_Summary.zip"
# zipfile_bed = "https://www.research-collection.ethz.ch/bitstream/handle/20.500.11850/434697/07_GlacierBed_SGI.zip"
# file_bed = "07_GlacierBed_SGI/B43-03_GlacierBed.tif"
#
# # download
# destdir = mkpath(joinpath(R20.datadir, "topography"))
# zipfile_bed = R20.download_file_to_dir(zipfile_bed, destdir)
# zipfile_outline = R20.download_file_to_dir(zipfile_outline, destdir)
# files_surf = [R20.download_file_to_dir(f, destdir) for f in files_surf]
# zipfile_CH = R20.download_file_to_dir(zipfile_CH, destdir)
# # merge tifs into one
# file_surf = joinpath(destdir, "surface_DEM.tif")
# run(`gdal_merge.py -of GTiff -o $file_surf $files_surf`)
#
# # unzip (requires unzip)
# run(`unzip -o $zipfile_outline -d $destdir`)
# shpfile_outline = joinpath(destdir, "01_Summary/Summary.shp")
#
# run(`unzip -o $zipfile_bed $file_bed -d $destdir`)
# file_bed = joinpath(destdir, file_bed)
#
# run(`unzip -o $zipfile_CH $folder_CH -d $destdir`) # why is this still in the folders?
# shpfile_CH = joinpath(destdir, "BOUNDARIES_2021/DATEN/swissBOUNDARIES3D/SHAPEFILE_LV95_LN02/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp")

# # download shapefiles of surface streams from glazio
# if !isdir(destdir * "/streams")
#     mkdir(destdir * "/streams")
#     for f in readdir(R20.glaziodir * "map_data/stream_coordinates")
#         R20.download_file("file://" * R20.glaziodir * "map_data/stream_coordinates/" * f,
#                            joinpath(destdir, "streams/" * f))
#     end
# end

# download borehole coordinates from glazio
# R20.download_file("file://" * R20.glaziodir * "map_data/BH_coordinates.csv",
#                       joinpath(destdir, "BH_coordinates.csv"))


# ------------------------------------------ #
#          Load data from data folder        #
# ------------------------------------------ #


fielddir = "../../data/geodata/field_data/"
ethdir = "../../data/geodata/research_collection_ethz/"
swisstopodir = "../../data/geodata/swisstopo/"

# shapefile glacier outline
shpfile_outline = ethdir * "outline_rhonegletscher/Summary.shp"
table = Shapefile.Table(shpfile_outline)
outline = nothing
for row in table
    if row.pk_sgi=="B43/03"
        global outline = Shapefile.shape(row)
    end
end

# shapefile CH boundary
shpfile_CH = swisstopodir * "outline_CH/swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp"
table = Shapefile.Table(shpfile_CH)
boundary = nothing
for row in table
    if row.NAME=="Schweiz"
        global boundary = Shapefile.shape(row)
    end
end

# shapefile streams
shpfile_stream1 = fielddir * "streams/stream_AM13_15.shp"
table = Shapefile.Table(shpfile_stream1)
stream1 = nothing
for row in table
    global stream1 = Shapefile.shape(row)
end

shpfile_stream2 = fielddir * "streams/stream_AM14.shp"
table = Shapefile.Table(shpfile_stream2)
stream2 = nothing
for row in table
    global stream2 = Shapefile.shape(row)
end

# bed dem
file_bed = ethdir * "B43-03_GlacierBed.tif"
bed = GDALarray(file_bed)[:,:,1] # the indexing brings it into memory

# surface dem
file_surf = swisstopodir * "surface_DEM.tif"
surf = GDALarray(file_surf)[1:5:end,1:5:end,1]
# thick = surf.-bed # does not work
# smooth it
function boxcar(A, window)
    out = similar(A)
    R = CartesianIndices(size(A))
    I1, Iend = first(R), last(R)
    for I in R
        n, s = 0, zero(eltype(out))
        I_ul = CartesianIndex(I1.I .* window)
        for J in CartesianIndices(UnitRange.(max(I1, I-I_ul).I , min(Iend, I+I_ul).I) )
            s += A[J]
            n += 1
        end
        out[I] = s/n
    end
    out
end
surf = boxcar(surf, 3)
#thick = boxcar(thick, 10)

# borehole coordinates
BH_file = fielddir * "BH_coordinates.csv"
BH_data = readdlm(BH_file, ','; skipstart=1)
BH_easting = BH_data[:,2] .*1e-3 # in km
BH_northing = BH_data[:,3] .*1e-3 # in km


#######
## Plotting
#
# https://discourse.julialang.org/t/labels-for-levels-in-contour-plot-plots-jl

# cl = :grey
# p = contour(bed, levels=2000:50:3000,  ls=:dash,
#             contour_labels=true,
#             #color = cgrad([:black,:black,:black]),
#             color = cgrad([cl, cl, cl]),
#             aspect_ratio = :equal,
#             colorbar=false
#             )
# contour!(surf, levels=2000:50:3000,
#          contour_labels=true,
#          color = cgrad([:black,:black,:black]),
#          colorbar=false,
#          )
# plot!(outline, fillalpha=0,
#       lw=4,
#       #color=:blue,
#       xlim=(2671454, 2673415),
#       ylim=(1160739, 1162721),
#       ylabel = "Northing (m)",
#       xlabel = "Easting (m)",
#       title="",
#       axes=:box,
#       )

# using PyPlot instead
# https://discourse.julialang.org/t/labels-for-levels-in-contour-plot-plots-jl/8266/3

cl = :grey
bed[bed.==-9999] .= NaN
levels_surface = 2000.0:20:3000
level_labels_surf = map(x -> "$(round(Int,x))", levels_surface)
level_ind_surf = 1:1:length(level_labels_surf)
levels_bed = 2000.0:50:3000
level_labels_bed = map(x -> "$(round(Int,x))", levels_bed)
level_ind_bed = 1:1:length(level_labels_bed)

Plt.figure(figsize=(14, 6))

# font sizes
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.size"] = 14 # general
fs_panellabels = 19 # for panel labels a-c

# positions and sizes of figures
xmin_CH, xmax_CH, ymin_CH, ymax_CH = 2468.0, 2851.0, 1064.0, 1307.0
xmin_Rhone, xmax_Rhone, ymin_Rhone, ymax_Rhone = 2671.2, 2675.6, 1158.8, 1167.8
ratio_CH = (ymax_CH-ymin_CH)/(xmax_CH-xmin_CH)
ratio_Rhone = (ymax_Rhone-ymin_Rhone)/(xmax_Rhone-xmin_Rhone)

left = 0.05
lpad = 0.082
bottom = 0.1
width = 0.12
heightCH = width*ratio_CH*14/6
heightRhone = width*ratio_Rhone*14/6

# surface topography
ax = Plt.axes([left, bottom, 0.9, 0.85])
cs = Plt.contour(surf.dims[1].val./1e3, surf.dims[2].val./1e3, surf',
                levels=levels_surface, colors="k", linewidths=1)
Plt.xlim((2671.5, 2673.5))
Plt.ylim((1160.8, 1162.8))

Plt.clabel(cs, levels_surface[level_ind_surf], fmt = Dict(zip(levels_surface[level_ind_surf], level_labels_surf[level_ind_surf])),
           inline =true,
           inline_spacing=1)

# bed topography
cs = Plt.contour(bed.dims[1].val./1e3, bed.dims[2].val./1e3, bed',
                levels=levels_bed, colors="grey", linewidths=1)
Plt.xlim((2671.5, 2673.5))
Plt.ylim((1160.8, 1162.8))
Plt.clabel(cs, levels_bed[level_ind_bed], fmt = Dict(zip(levels_bed[level_ind_bed], level_labels_bed[level_ind_bed])),
           inline = 1,
           inline_spacing=1)

# glacier outline
outl = hcat([[p.x,p.y] for p in outline.points]...);
Plt.plot(outl[1,:]/1e3, outl[2,:]/1e3, lw=3)

# BH/AM coordinates
p1 = Plt.scatter(BH_easting[[1,2,4]], BH_northing[[1,2,4]], 30, "k", "o") # unactivated BH
p2 = Plt.scatter(BH_easting[[3,5]], BH_northing[[3,5]], 50, "k", "x") # AM

Plt.annotate("BH11", (BH_easting[1]-0.085, BH_northing[1]-0.005))
Plt.annotate("BH12", (BH_easting[2]-0.085, BH_northing[2]-0.01))
Plt.annotate("AM13/BH13", (BH_easting[3]-0.15, BH_northing[3]-0.005))
Plt.annotate("BH14", (BH_easting[4]+0.02, BH_northing[4]-0.015))
Plt.annotate("AM15/BH15", (BH_easting[5]-0.15, BH_northing[5]))

# stream
str1 = hcat([[p.x,p.y] for p in stream1.points]...);
str2 = hcat([[p.x,p.y] for p in stream2.points]...);
p4 = Plt.plot(str1[1,:]/1e3, str1[2,:]/1e3, "k--", lw=1)
Plt.plot(str2[1,:]/1e3, str2[2,:]/1e3, "k--", lw=1)
p3 = Plt.scatter(str1[1,1]/1e3, str1[2,1]/1e3, 30, "k", "<")
Plt.scatter(str2[1,1]/1e3, str2[2,1]/1e3, 30, "k", "<")

# workaround to make legend for 2D lines
L2D = Plt.matplotlib.lines.Line2D
custom_lines = [L2D([0], [0], color="k", linewidth=1, linestyle="--"),
                L2D([0], [0], color="k", linewidth=1),
                L2D([0], [0], color="grey", linewidth=1),
                L2D([0], [0], linewidth=3)]
Plt.legend((custom_lines[1], custom_lines[2], custom_lines[3], custom_lines[4],
            p1, p2, p3),
            ("Surface streams", L"\mathrm{Surface\,\,contour\,\,lines,\,\,20\,m\,\,spacing}", L"\mathrm{Bed\,\,contour\,\,lines,\,\,50\,m\,\,spacing}", "Glacier outline",
            "Not activated boreholes", "Artificial moulins", "Natural moulins"),
            framealpha = 1.0,
            loc = "lower right",
            ncol = 2,
            fontsize = 15,
            )

ymin, ymax = 1160.63, 1161.23
xmin, xmax = 2671.55, 2672.8
Plt.axis("equal")
Plt.axis("square")
Plt.xlim((xmin, xmax))
Plt.ylim((ymin, ymax))
Plt.xlabel("Easting (km)")
Plt.ylabel("Northing (km)")
# Plt.tick_params(labelleft=false, left=false, labelbottom=false, bottom=false)
Plt.text(0.88, 0.98, L"\bf{c}", transform=ax.transAxes, ha="left", va="top", fontsize=fs_panellabels+2) #
# Plt.hlines(ymax-0.25*height1/(height1+height2)*(ymax-ymin), xmin+0.895, xmin+0.995,
#           "k", linewidth=4)
# Plt.annotate(L"\mathrm{100\,m}", (xmin+0.92, ymax-0.22*height1/(height1+height2)*(ymax-ymin)), fontsize=12)
# Plt.arrow(xmin+0.945, ymin+0.38, 0, 0.055,
#          facecolor="black",
#          width=0.002,
#          head_width=0.01,
#          head_length=0.015
#          )
# Plt.annotate("N", (xmin+0.925, ymin+0.39), fontsize=28)
Plt.minorticks_off()


# glacier outline alone
ax = Plt.axes([left+lpad, bottom+(0.85-heightCH-heightRhone)*1/3, width, heightRhone]) #figsize=(12, 7))
Plt.plot(outl[1,1:1364]/1e3, outl[2,1:1364]/1e3, lw=2, zorder=1) # exclude some points at the end which are within the glacier and make outline ugly
Plt.hlines([ymin, ymax], xmin, xmax, "k", lw=2, zorder=2)
Plt.vlines([xmin, xmax], ymin, ymax, "k", lw=2, zorder=2)
Plt.axis("equal")
Plt.axis("square")
Plt.tick_params(labelleft=false, left=false, labelbottom=false, bottom=false)
Plt.xlim([xmin_Rhone, xmax_Rhone])
Plt.ylim([ymin_Rhone, ymax_Rhone])
Plt.text(xmin+0.5*(xmax-xmin), ymax+0.3*(ymax-ymin), "c", fontsize=fs_panellabels, ha="center")
Plt.text(0.06, 0.98, L"$\bf{b}$", transform=ax.transAxes, ha="left", va="top", fontsize=fs_panellabels)
Plt.hlines(ymin_Rhone+0.15*(ymax_Rhone-ymin_Rhone), xmin_Rhone+0.6*(xmax_Rhone-xmin_Rhone), xmin_Rhone+0.6*(xmax_Rhone-xmin_Rhone)+1.0,
           "k", linewidth=4)
Plt.annotate(L"\mathrm{1\,km}", (xmin_Rhone+0.6*(xmax_Rhone-xmin_Rhone), ymin_Rhone+0.18*(ymax_Rhone-ymin_Rhone)), fontsize=12)
Plt.arrow(0.55, 0.8, 0., -0.1, width=0.01, head_width=0.05, head_length=0.03, transform=ax.transAxes, facecolor="black")
Plt.text(0.65, 0.75, "ice flow", transform=ax.transAxes, ha="left", va="center", rotation="vertical")

# CH outline alone
ax = Plt.axes([left+lpad, 0.95-(0.85-heightCH-heightRhone)*2/3-heightCH, width, heightCH])
bCH = hcat([[p.x,p.y] for p in boundary.points]...);
Plt.plot(bCH[1,1:end-414]/1e3, bCH[2,1:end-414]/1e3, "k") # exclude some points, same reason as in glacier outline
Plt.scatter(BH_easting[1], BH_northing[1])
Plt.axis("equal")
Plt.tick_params(labelleft=false, left=false, labelbottom=false, bottom=false)
Plt.xlim([xmin_CH, xmax_CH])
Plt.ylim([ymin_CH, ymax_CH])
ymin_CH, ymax_CH = Plt.ylim()
xmin_CH, xmax_CH = Plt.xlim()
Plt.text(0.06, 1-0.02*heightRhone/heightCH, L"$\bf{a}$", transform=ax.transAxes, ha="left", va="top", fontsize=fs_panellabels)
Plt.text(BH_easting[1]-15, BH_northing[1]+30.0, "b", fontsize=fs_panellabels)

Plt.gcf()
