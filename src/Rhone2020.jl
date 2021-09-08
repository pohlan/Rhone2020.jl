__precompile__(false)
module Rhone2020


# ----------------------------------------------------------------- #
#                                                                   #
#                              Packages                             #
#                                                                   #
# ----------------------------------------------------------------- #

using DelimitedFiles
using Dates
using StringEncodings
using MonteCarloMeasurements
using LsqFit
using Statistics
using PyPlot


# ----------------------------------------------------------------- #
#                                                                   #
#                        Define constants                           #
#                                                                   #
# ----------------------------------------------------------------- #

# physical constants
const g = 9.81 # gravitational acceleration in m/s^2
const rhow = 1000.0 # density of water in kg/m^3
const rhoi = 910.0 # density of ice, kg/m^3
const L = 334e3 # latent heat of fusion for water, J/kg
const ct_pure = -7.4e-8 # pressure melt coefficient for pure water, K/Pa
const ct_air = -9.8e-8 # pressure melt coefficient for air saturated water, K/Pa
const cw = 4217.7 # heat capacity of water, J/(K*kg)
const mu = 1.787e-3 # dynamic viscosity of water in Pa s
const A = 2.4e-24 # creep parameter, 1/(s*Pa^3), for 0°C
const n = 3 # flow law exponent
const atm = 101.325e3 # standard pressure in Pa (where melting temperature is 0°C)

# directory where the data resides
const datadir = joinpath(@__DIR__, "../../data/")

# export to make easily accessible outside of module Rhone2020.jl
export rhow, g, mu, cw, ct, ct_pure, ct_air, n_partcl, set_n_partcl

"Function to (re)set number of particles used in MC uncertainty propagation"
function set_n_partcl(n)
    # Sample size for uncertainty propagation with MonteCarloMeasurements.jl
    global n_partcl = n # 20000 used for paper because Fig. 4 looks nicer then (but takes some time), 2000 for plots in data collection
    # define ct as a uniform distribution between the two literature values
    global ct = Particles(n_partcl, Uniform(ct_air, ct_pure))
    return nothing
end
set_n_partcl(1000)

# set resolution for all figures
pixel_resolution = 300.      # dpi = dots per inch

# ----------------------------------------------------------------- #
#                                                                   #
#                Functions to download and read data                #
#                                                                   #
# ----------------------------------------------------------------- #

"""
    download_file(url, destination_file, force_download=false)
Download a file, if it has not been downloaded already.
Input:
- url -- url for download
- destination_file -- path of the file to store the download
- force_download -- force the download, even if file is present
Output: nothing
"""
function download_file(url, destination_file; force_download=false)
    if !isfile(destination_file) || force_download
        download(url, destination_file)
    end
    return nothing
end

"""
Download to the directory and figure out filename automatically.
"""
function download_file_to_dir(url, destination_dir=datadir; force_download=false)
    destination_file = joinpath(destination_dir, splitdir(url)[2])
    if !isfile(destination_file) || force_download
        download(url, destination_file)
    end
    return destination_file
end



"""
    read_WTW(filename)
This function reads a file from the WTW conductivity sensor and
returns:
Dict with keys: :t [date-time stamp], :cond [μS/cm], :temp [C]
Note, that the input file usually contains several traces.  Split them up with
`split_conductivity_data`.
"""
function read_WTW(filename)
    if !isfile(filename)
        error("Filename $filename is not a file!")
    end
    if !startswith(splitdir(filename)[2], "AD")
        warn("Read in a file starting with `AD`.  (The `AC` files use a comma for the decimal point.")
    end
    d, head = readdlm(filename, ';', header=true)
    out = Dict{Symbol,Any}() # key has the be Symbol, value can be anything
    # time 12.08.2016 13:36:58
    fmt = "d.m.y H:M:S"
    out[:t] = [DateTime(dd, fmt) for dd in d[:,4]]
    # conductivity
    out[:cond] = d[:,5]
    units = d[:,6]
    @assert all(units.=="\xb5S/cm") "Units not in μS/cm!"
    # temp
    out[:temp] = d[:,8]
    # ID
    out[:id] = d[:,3]

    # purge any records which are simultaneous (not sure why this happens with the WTW)
    purge = findall(diff(out[:t]).==Second(0))
    for p in reverse(purge)
        deleteat!(out[:t], p)
        deleteat!(out[:cond],p)
        deleteat!(out[:temp], p)
        deleteat!(out[:id], p)
    end

    return out
end

"""
         read_Keller(filename;
                     presshead="P1",
                     condhead="ConRaw",
                     temphead="TOB1",
                     skipstart=8,
                     )
Reads a Keller pressure/CTD sensor.  However, you probably want to use
- `read_Keller_DCX22_CTD`,
- `read_Keller_DCX22` and
- `read_Keller_DC22`
Returns a dict with keys as appropriate:
- :t [date-time stamp]
- :cond [μS/cm]
- :temp [C]
- :press [m H2O]
"""
function read_Keller(filename;
                     presshead="P1",
                     presshead2="P2",
                     condhead="ConRaw",
                     temphead="TOB1",
                     temphead_PT = "T",
                     temphead2="TOB2",
                     skipstart=8,
                     )
    d,h = readdlm(open(filename, enc"UTF-16","r"), '\t', skipstart=skipstart, header=true)
    h = h[:] # h is a 1x2 matrix, change to a vector

    out = Dict{Symbol,Any}()
    # find date-time rows
    id, it = findfirst(h.=="Date"), findfirst(h.=="Time")
    # time 12.08.2Q309[d]016 13:36:58
    fmtd, fmtt = "d/m/y", "H:M:S"
    out[:t] = [Date(dd, fmtd) + Time(tt, fmtt) for (dd,tt) in zip(d[:,id], d[:,it])]

    #
    for (head, key) in [(presshead, :press),
                        (presshead2, :airpress),
                        (condhead, :cond),
                        (temphead, :temp_TOB),
                        (temphead_PT, :temp_PT),
                        (temphead2, :airtemp)]
        i = findfirst(h.==head) # see if there is one
        tmp = Float64[]
        if i!==nothing
            out[key] = [s=="" ? missing :
                        s isa AbstractString ? parse(Float64, replace(s, ","=>".")) : Float64(s) for s in d[:,i]]
            # convert mS/cm to μS/cm
            if key==:cond
                out[:cond] = out[:cond].*1000
            end
            # convert pressure to Pa
            if key==:press
                out[:press] = out[:press].*1e2
            end
        end
    end

    # check lengths and remove all "missing"
    l = length(out[:t])
    topurge = []
    for v in values(out)
        @assert length(v)==l
        append!(topurge, findall(v.===missing))
    end
    topurge = sort(unique(topurge))
    for (k,v) in out
        deleteat!(v, topurge)
        if k!=:t
            out[k] = Float64.(v) # make the vector an
        end
    end
    return out
end

read_Keller_DCX22_CTD(filename) = read_Keller(filename)
read_Keller_DCX22_CTD_tc(filename) = read_Keller(filename, condhead="ConTc") # temp corrected
read_Keller_DCX22(filename) = read_Keller(filename)
read_Keller_DC22(filename) = read_Keller(filename,
                                         presshead="p1/Pa",
                                         temphead="t1/\xb0C",
                                         skipstart=24
                                         )
read_Keller_stage(filename) = read_Keller(filename)

"""
        load_DCX22_CTD():

Loads all data of DCX22_CTD sensors

Output: one dictionary for each ctd containing an entry for each day (e.g. ctd145["1408"] => data of ctd 205145 from 14.08.2020)
- 1st dictionary: ctd 205145-10mH2O
- 2nd dictionary: ctd 205309-100mH2O
- 3rd dictionary: ctd 207265-300mH2O
"""
function load_DCX22_CTD()
    ctd145 = Dict()
    ctd309 = Dict()
    ctd265 = Dict()
    for file in readdir(datadir * "DCX22_CTD/")
        strings = split(file, "_")
        key = join([strings[2], strings[3]])
        if startswith(file, "205145")
            ctd145[key] = read_Keller_DCX22_CTD(datadir * "DCX22_CTD/" * file)
        elseif  startswith(file, "205309")
            ctd309[key] = read_Keller_DCX22_CTD(datadir * "DCX22_CTD/" * file)
        elseif startswith(file, "207265")
            ctd265[key] = read_Keller_DCX22_CTD(datadir * "DCX22_CTD/" * file)
        end
    end
    ctd145["number"]="145"
    ctd309["number"]="309"
    ctd265["number"]="265"
    return ctd145, ctd309, ctd265
end

"""
        load_DCX22_pressure():

Loads all data of DCX22_pressure sensor

Output: one dictionary containing an entry for each file
"""
function load_DCX22_pressure()
    pressure = Dict()
    for file in readdir(datadir * "DCX22_pressure/")
        strings = split(file, "_")
        key = join([strings[2], strings[3]])
        pressure[key] = read_Keller_DCX22(datadir * "DCX22_pressure/" * file)
    end
    return pressure
end

"""
        load_stage_pressure():

Loads all data of stage pressure sensor

Output: dictionary with one time series for each data set
"""
function load_stage_pressure()
    pressure = Dict()
    for file in readdir(datadir * "Stage_pressure/")
        strings = split(file, "_")
        key = join([strings[2], strings[3]])
        current = read_Keller_stage(datadir * "Stage_pressure/" * file)
        if haskey(pressure,:t)
            for k in keys(current)
                pressure[k]=append!(pressure[k],current[k])
            end
        else
            pressure = current
        end
    end
    return pressure
end


# ----------------------------------------------------------------- #
#                                                                   #
#                  Functions to derive Q, S, f, etc.                #
#                                                                   #
# ----------------------------------------------------------------- #

"""
Converts ml added to bucket to a concentration (g/l == kg/m^3).
Input:
- ml -- how many mililiters were added
- solution -- the concentration of the calibration solution (kg/m^3 == g/l)
- bucketsize -- the size of the bucket/bottle to which the solution was added (l)
Output:
- concentration (kg/m^3 == g/l)
"""
function ml_to_concentration(ml, solution, bucketsize)
    mass = ml/1e3 * solution # salt mass added to bucket (g)
    return mass/bucketsize # concentration in g/l (== kg/m^3)
end

"""
    curve_fit_MCMeasurements(model, x, y, p0)

Does LsqFit.curve_fit with MonteCarloMeasurements Vectors.  Same input as `curve_fit`.
"""
function curve_fit_MCMeasurements(model, x, y, p0)
    if eltype(x)<:Particles && eltype(y)<:Particles
        len = length(x[1].particles)
        pout = [eltype(p0)[] for p in p0]
        for i=1:len
            p = curve_fit(model, [xx.particles[i] for xx in x], [yy.particles[i] for yy in y], p0)
            @assert p.converged
            for j=1:length(p0)
                push!(pout[j], p.param[j])
            end
        end
    elseif eltype(y)<:Particles
        len = length(y[1].particles)
        pout = [eltype(p0)[] for p in p0]
        for i=1:len
            p = curve_fit(model, x, [yy.particles[i] for yy in y], p0)
            @assert p.converged "$i"
            for j=1:length(p0)
                push!(pout[j], p.param[j])
            end
        end
    elseif eltype(x)<:Particles
        len = length(x[1].particles)
        pout = [eltype(p0)[] for p in p0]
        for i=1:len
            p = curve_fit(model, [xx.particles[i] for xx in x], y, p0)
            @assert p.converged
            for j=1:length(p0)
                push!(pout[j], p.param[j])
            end
        end
    else
        error()
    end
    return [Particles(po) for po in pout]
end

"""
    fit_calibration(bucketsize, solution, calis...)
Fits a line of best fit through the calibration data going through the origin.
Returns the function of this line: `f(cond-cond_at_0) -> concentration`.
Also prints the parameters values +/- 95% confidence intervals.
Uses the package LsqFit: https://github.com/JuliaOpt/LsqFit.jl
"""
function fit_calibration(bucketsize, solution, error_cond, calis...)
    # subtract readout at zero concentration
    for c in calis
        @assert c[1,1]==0 "First row needs to be zero reading!"
        c[:,2] .= c[:,2].-c[1,2]
    end
    # concatenate all calibrations
    cali = vcat(calis...)
    conc = ml_to_concentration(cali[:,1], solution, bucketsize)
    delta_readout = Particles.(n_partcl, Normal.(cali[:,2], error_cond))
    # Fit line using https://github.com/JuliaOpt/LsqFit.jl
    fn(delta_readout, p) = p[1]*delta_readout # p[1]==a, p[2]==b
    para_weights = [0.5] # equal weights to parameters

    fit = curve_fit_MCMeasurements(fn, delta_readout, conc, para_weights)
    #errors = margin_error(fit, 1-0.95)
    #println("""
    #Estimated linear fit: f(delta_cond) = a*conc with
    # a = $(round(fit.param[1],sigdigits=3))±$(round(errors[1],sigdigits=3))
    #""")
    return (delta_readout) -> fn(delta_readout, fit)
end;

"""
    boxcar(A::AbstractArray, window, [, weights, keepmask])
    boxcar(A::AbstractArray, windows::Tuple, [, weights, keepmask])
    boxcar(A::AbstractArray, window::AbstractArray, [, weights, keepmask])
Boxcar filter.  The two argument call skips NaNs.  The three & four
argument call uses weights and propagates NaNs, it can be a lot faster.
Smoothing occurs over +/-window indices, 0 corresponds to no smoothing.
The window can be specified as:
- integer, for a symmetric window in all dimensions
- a tuple to give lower and upper windows
- a tuple of tuples to give different lower and upper windows for all dimensions
- a array of size(A) for a different, symmetric window at each point.
Weights, if given, will use those relative weights for averaging.  Note that
points which have value==NaN and weight==0 will not poison the result.
No average is calculated for points where keepmask==true, instead
their original value will be kept.
Notes:
- For the weights it may be faster to use non-Bool arrays nor BitArrays,
  say Int8.  Note that NaNs where weight==0 will not poison the result.
- Also works for Vectors.
From http://julialang.org/blog/2016/02/iteration
"""
boxcar(A::AbstractArray, window) = boxcar(A, (window,window))
function boxcar(A::AbstractArray, windows::Tuple)
    window_lower, window_upper = windows
    out = similar(A)
    R = CartesianIndices(size(A))
    I1, Iend = first(R), last(R)
    I_l = CartesianIndex(I1.I.*window_lower)
    I_u = CartesianIndex(I1.I.*window_upper)
    for I in R # @inbounds does not help
        out[I] = NaN
        n, s = 0, zero(eltype(out))
        for J in CartesianIndices(UnitRange.(max(I1, I-I_l).I , min(Iend, I+I_u).I) ) # used to be
        # for J in max(I1, I-I_l):min(Iend, I+I_u) ## in Julia 1.1
            if !isnan(A[J])
                s += A[J]
                n += 1
            end
        end
        out[I] = s/n
    end
    out
end

"""
    discharge(trace,date,indices,mass,error_mass,error_cond)

Calculates the discharge: Q=mass/(integral of concentration)

Input:
- trace -- ctd file, e.g. ctd309_2008
- date -- string in format ddmm, e.g. "0808"
- indices -- array with start and end points of salt tracer signal
- mass -- array with masses of each injection e.g. m = [0.1,0.1,0.1,0.1]
- error_mass -- error of mass measurements in kg
- error_cond -- error of conductivity measurements in μS/cm

Output:
- vector of discharges in m^3/s
"""
function discharge(trace, date, indices, mass, error_mass, error_cond)
    Qs = []
    t, cond = trace[date][:t], trace[date][:cond] # time and conductivity data of the whole day
    for k in 1:size(indices[1],2) # loop through injections of one day
        i1 = indices[1][k]
        i2 = indices[2][k]
        if date == "2008" || date == "2108"
            delta_cond = Particles.(n_partcl, Normal.(
                         cond[i1:i2] .- mean(cond[[i1;i2]]), # readings do not return to background value for a long time
                         error_cond))
        else
            delta_cond = Particles.(n_partcl, Normal.(
                         cond[i1:i2] .- mean(cond[[i1-30:i1;i2:i2+30]]),
                         error_cond))
        end
        conc = trace["cond2conc"].(delta_cond)
        dt = Dates.Second(t[2]-t[1])/Dates.Second(Time(0,0,1)) # sampling interval in seconds
        integral_conc = sum(conc*dt)
        Q = Particles.(n_partcl, Normal.(
            mass[date][k], error_mass)
            )  /integral_conc # discharge in m^3/s
        push!(Qs, Q)
    end
    return Qs
end

"""
    flow_speed(date, indices, t0, depths, trace1, trace2, error_dz, error_dt)

Calculates the average flow speed between two sensors and between the upper sensor and the surface.

Input:
- date -- string in format ddmm, e.g. "0808"
- indices -- tuple of 4 arrays with start and end points of salt tracer signal in each of the ctds
- t0 -- dictionary with indexes of injection times for each date, e.g. t0["0808"] = ([16806 18186 19606 21906])
- depths -- array with depths of the ctds in m; first the upper, second the lower ctd depth
- ctd1 -- dictionary of upper ctd holding conductivity measurements and time
- ctd2 -- dictionary of lower ctd holding conductivity measurements and time
- error_dz -- estimated error of depths difference in m
- error_dt -- error of time difference between manually determined injection time and peak of conductivity, in s
    (error for time difference between two ctds is fixed as half the sampling interval)

Output:
- tinj -- vector with the injection times
- V -- vector with average flow velocities between sensors, one value for each injection, in m/s
- V_surf -- vector with average flow velocities between surface and upper sensor, one value for each injection, in m/s
"""
function flow_speed(ctd1, ctd2, date, indices, t0, error_dt)
    tinj = DateTime[]
    Vbot_peak = []
    Vtop_peak = []
    Vbot_HM = []
    Vtop_HM = []

    for k in 1:size(indices[1],2) # loop through injections of one day
        i1 = indices[1][k]
        i2 = indices[2][k]
        i3 = indices[3][k]
        i4 = indices[4][k]

        cond1, cond2 = ctd1[date][:cond], ctd2[date][:cond] # conductivity vectors of the whole day
        d1, d2 = ctd1[date][:depth], ctd2[date][:depth] # sensor depths, m

        t = ctd1[date][:t] # time vector of the whole day (same for both ctds)
        dt = Dates.Second(t[2]-t[1])/Dates.Second(Time(0,0,1)) # sampling interval in seconds

        # compute velocity from time difference between peaks
        ip_inj = t0[date][k] # index of injection
        ip1 = findmax(cond1[i1:i2])[2]+i1-1 # index of peak in upper ctd
        ip2 = findmax(cond2[i3:i4])[2]+i3-1 # index of peak in lower ctd
        vbot_peak = (d2-d1) ./                                             # velocity between upper and lower ctd in m/s
                    Particles.(n_partcl, Normal.((ip2-ip1)*dt, 0.5*dt))
        vtop_peak = d1 ./                                                  # velocity between surface and upper ctd in m/s
                    Particles.(n_partcl, Normal.((ip1-ip_inj)*dt, error_dt))

        push!(Vbot_peak, vbot_peak)
        push!(Vtop_peak, vtop_peak)

        # compute velocity from harmonic mean times
        integral1 = sum(cond1[i1:i2] * dt) # integral of conductivity
        integral2 = sum(cond2[i3:i4] * dt)
        pdf1 = cond1[i1:i2] / integral1 # normalised distribution, corresponds to probability density function
        pdf2 = cond2[i3:i4] / integral2
        idx1 = i1:i2
        idx2 = i3:i4
        idx1_HM = (sum(pdf1 ./ idx1 * dt))^(-1) # harmonic mean index
        idx2_HM = (sum(pdf2 ./ idx2 * dt))^(-1)
        vbot_HM = (d2-d1) ./                                                            # velocity between upper and lower ctd in m/s
                            Particles.(n_partcl, Normal.((idx2_HM-idx1_HM)*dt, 0.5*dt))
        vtop_HM = d1 ./                                                                 # velocity between surface and upper ctd in m/s
                            Particles.(n_partcl, Normal.((idx1_HM - ip_inj) * dt, error_dt))

        push!(Vtop_HM, vtop_HM)
        push!(Vbot_HM, vbot_HM)
        push!(tinj,t[Int(floor(mean([idx1_HM, idx2_HM])))])

    end
    return tinj, Vtop_peak, Vbot_peak, Vtop_HM, Vbot_HM
end

"""
    hydraulic_gradient(ctd1, ctd2, mid, date, indices, error_p)

Calculates the hydraulic gradient dphi/dz = dp/dz - rhow * g

Input:
- ctd1 -- dictionary of upper ctd holding pressure and time
- ctd2 -- dictionary of lower ctd holding pressure and time
- mid -- dictionary holding derived values of test-section between ctd1 and ctd2, in particular the length of the test-section, dz
- date -- the day for which the hydraulic gradient is be calculated, e.g. "0908" for 09 August

Output:
- dphi_dz -- hydraulic gradient between the two sensors at time points of the experiments, Pa/m
- dp -- pressure difference between the two sensors at time points of the experiments, Pa
- p_mean -- mean pressure between the two sensors at time points of the experiments, Pa
"""
function hydraulic_gradient(ctd1, ctd2, mid, date, indices, error_p)
    ctdup, ctdlow = ctd1[date], ctd2[date]
    dz = mid[date][:dz]  # vertical distance between ctds in m

    if ctd1["number"] == "145"
        error_p_up = 10.0 * error_p
        error_p_low = 100.0 * error_p
    else
        error_p_up = 100.0 * error_p
        error_p_low = 300.0 * error_p
    end

    dp = []
    p_mean = []
    for k in 1:length(indices[1])
        ix = indices[1][k]:indices[4][k] # start salt signal upper ctd : end salt signal lower ctd
        dpress = Particles.(n_partcl, Normal.(ctdlow[:press][ix], error_p_low)) .- # hydraulic gradient in kPa
                 Particles.(n_partcl, Normal.(ctdup[:press][ix], error_p_up))
        idx_stage = indexin(DateTime.(ctdup[:t][ix]), DateTime.(ctdup[:t_dpress])) # stage sensor (with air pressure) has larger sampling interval
        if sum(idx_stage .!= nothing) == 0
            k = 1
            while sum(idx_stage .!= nothing) <= 1
                idx_stage = indexin(DateTime.(ctdup[:t][ix[1]-k:ix[end]+k]), DateTime.(ctdup[:t_dpress]))
                k = k+1
            end
        end
        idx_stage = idx_stage[idx_stage.!==nothing]
        push!(dp, only(mean(dpress, dims=1))) # average pressure difference over integration bounds
        push!(p_mean, only(mean([only(mean(Particles.(n_partcl, Normal.(ctdlow[:dpress][idx_stage], error_p_low)), dims=1)),
                            only(mean(Particles.(n_partcl, Normal.(ctdup[:dpress][idx_stage], error_p_up)), dims=1))],
                            dims=1)))
    end

    # hydraulic gradient
    dphi_dz = dp/dz .- rhow*g

    return dphi_dz, dp, p_mean

end

"""
    friction(mid, date, vtype)

calculates the Darcy-Weisbach an manning friction

Input:
- mid -- dictionary holding derived values of test-section between ctd1 and ctd2, in particular dphi/dz, Q, S, type of v calculation (based on peaks vs. harmonic mean)
- date -- day where frictions are to be calculated,

"""
function friction(mid, date, vtype)
    dphi_dz = mid[date][:dphi_dz]
    Q, v, S = mid[date][:Q], mid[date][vtype], mid[date][:S]

    # friction factor, unitless
    f = 4*abs.(dphi_dz).*sqrt.(Q/pi)./(rhow*v.^(2.5))

    # manning roughness, s/m^{1/3}
    RH = 0.5*sqrt.(S./pi) # hydraulic radius, m
    n_manning = sqrt.(RH.^(1/3).*f./(8*g))

    return f, n_manning
end

"""
Calculate Reynolds number: Re = rho_w * v * D / mu, where
rho_w: density of water, v: flow speed, D: diameter, mu: dynamic viscosity of water

"""
function Reynolds_number(mid, date, vtype)
    S, v = mid[date][:S], mid[date][vtype]
    Re = 2*rhow/ mu .* v .* sqrt.(S./pi)
    return Re
end


# ----------------------------------------------------------------- #
#                                                                   #
#              Functions used for size evolution models             #
#                                                                   #
# ----------------------------------------------------------------- #

"""
    closure_rate(ctd1, ctd2, mid, date)

Calculates closure rates according to viscous deformation of ice (Glen's flow law)

Input:
- ctd1 -- dictionary of upper ctd holding depth
- ctd2 -- dictionary of lower ctd holding depth
- mid -- dictionary holding derived values of test-section between ctd1 and ctd2, in particular S and p_mean
- date -- the day for which the opening/closure rates are to be calculated, e.g. "0908" for 09 August

Output:
- closure -- closure rate, m^2/s
"""
function closure_rate(ctd1, ctd2, mid, date)
    # data
    S = mid[date][:S] # cross-sectional area, m^2
    p_mean = mid[date][:p_mean] # mean pressure, Pa
    d1, d2, = ctd1[date][:depth], ctd2[date][:depth]

    # closure rate due to ice creep
    N = rhoi * g * mean([d1, d2], dims=1) .- p_mean # effective pressure = p_ice - p_water, Pa
    closure = 2.0/(n^n) * A * S .* abs.(N).^(n-1) .* N # m^2/s

    return closure
end

"""
    model_S(mid, date, dT_dz, S0)

Models the size evolution of S for given temperature gradients.

Inputs:
- mid -- dictionary holding derived values of test-section between ctd1 and ctd2, in particular Q, dphi/dz and the injection times
- date -- the day for which the model is to be run, e.g. "0908" for 09 August
- dT_dz -- temperature gradient constant in time; can eighter be a Float64 (used in MCMC of free-gradient model) or Particles (used in ct-model)
- S0 -- initial cross-sectional area, at time of first experiment; can eighter be a Float64

Outputs:
- S_model -- modelled time evolution of the cross-sectional areas
- dS_dt -- total opening rates for each experiment
- sensible_heat -- parts of the opening rates that are due to sensible heat (without dissipation)
(The type of the outputs corresponds to the types of dT_dz and S0)
"""
function model_S(mid, date, dT_dz, S0)
    if typeof(dT_dz) == Float64 # for MCMC runs, only work with mean value (neglect uncertainties in Q and dphi/dz)
        Q, dphi_dz, t = pmean.(mid[date][:Q]), pmean.(mid[date][:dphi_dz]), mid[date][:t_inj]
    else # ct-model, use Particle distributions for propagation of uncertainties
        Q, dphi_dz, t = mid[date][:Q], mid[date][:dphi_dz], mid[date][:t_inj]
    end
    dt = (t - t[1])*1e-3/Dates.Millisecond(1) # time steps (time between two experiments)

    # opening rates
    dissipation = Q./(rhoi*L) .* abs.(dphi_dz)
    sensible_heat = Q./(rhoi*L) .* (-dT_dz*cw*rhow)
    dS_dt = dissipation + sensible_heat # m^2/s

    S_model = zeros(typeof(dS_dt[1]), length(dS_dt))
    S_model[1] = S0
    for i in 2:length(dS_dt)
        S_model[i] = only(S_model[i-1] .+ mean([dS_dt[i-1], dS_dt[i]], dims=1)*diff(dt)[i-1])
    end

    return S_model, dS_dt, sensible_heat
end

"""
    make_logpdf(date, mid, dTdzrange)

Returns a logarithmic probability density function as input for the free-gradient model

Inputs:
- date -- the day for which the model is to be run, e.g. "0908" for 09 August
- mid -- dictionary holding derived values of test-section between ctd1 and ctd2, in particular S, Q, dphi/dz and the injection times
- dTdzrange -- range of uniform distribution for prior in temperature gradient

Outputs:
- logarithmic pdf function
    Inputs: theta -- fitting parameters of the MCMC, dT/dz and initial S
    Outputs: L+prior -- logarithmic pdf
             (S_model, dS_dt, sensible_heat) -- outputs of model_S()
"""
function make_logpdf(date, mid, dTdzrange)
    dTdzmin, dTdzmax = dTdzrange
    S_data = pmean.(mid[date][:S])
    std_S = pstd.(mid[date][:S])

    return function(theta)
        dT_dz, S0 = theta

        if dTdzmin < dT_dz < dTdzmax # prior

            model_output = model_S(mid, date, dT_dz, S0)

            if typeof(model_output) == Array{Float64,1}
                S_model = pmean.(model_output)
            else
                S_model = pmean.(model_output[1])
            end
            # Gaussian likelihood function
            L = -0.5*sum((S_data.-S_model).^2 ./std_S.^2)
            # Gaussian prior probability for S0
            prior = -0.5*(S_data[1]-S0)^2/std_S[1]^2

            return L+prior, model_output # posterior, blob
        else
            return -Inf, nothing
        end
    end
end


# ----------------------------------------------------------------- #
#                                                                   #
#                        Plotting functions                         #
#                                                                   #
# ----------------------------------------------------------------- #

"""
Plot calibrations of salt concentrations vs. measured conductivity
"""
function plot_conc_cali(calis, ctds, error_m, error_cond)
    solution = 1
    bucketsize = 1

    # font size
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12

    figure(figsize=(12,3), dpi=pixel_resolution)
    for (c, number) in enumerate([309, 265, 145])
        # which ctd to take
        ctd = ctds[findfirst(in.(("number" => string(number)), ctds))]
        if number == 145
            days = ["23.06.","20.08."] # less calibrations with CTD-145
        else
            days = ["05.08.", "10.08.", "13.08.", "20.08."]
        end

        subplot(1, length(calis), c)

        # scatter plot of calibration readouts
        for n in 1:length(calis[number]) # loop through days where calibrations were conducted
            conc = ml_to_concentration(Particles.(n_partcl, Normal.(calis[number][n][:,1], error_m)),
                                       solution, bucketsize)
            cond = Particles.(n_partcl, Normal.(calis[number][n][:,2], error_cond)) .-
                   Particles.(n_partcl, Normal.(calis[number][n][1,2], error_cond))
            errorbar(pmean.(conc), pmean.(cond), xerr=pstd.(conc), yerr=pstd.(cond), label=days[n], fmt="+")
        end

        # plot the line of best fit:
        func = ctd["cond2conc"]
        cond_bestfit = range(0,stop=maximum(maximum.(calis[number])),length=100)
        plot(pmean.(func.(cond_bestfit)), cond_bestfit, color="black", label="best fit")
        xlabel("salt concentration [g/l]")

        title("CTD-" * ctd["number"])
        legend(handlelength=1.0)

        if c == 1
            ylabel("sensor conductivity [μS/cm]")
        end
    end
    gcf()
end

"""
Plot zero degree calibration of temperature sensors
"""
function plot_temp_cali(ctds, idx, date; dTOBs=[], dTs=[])
    # time axis formatting
    majorformatter = matplotlib.dates.DateFormatter("%H:%M");
    if date == "2306"
        majorlocator = matplotlib.dates.MinuteLocator(byminute=(10, 30, 50))
    end

    # draw figure
    figure(figsize=(12,3), dpi=pixel_resolution)
    for (k, (ctd, i, dTOB, dT)) in enumerate(zip(ctds, idx, dTOBs, dTs))
        subplot(1, length(ctds), k)
        plot(ctd[date][:t][i], ctd[date][:temp_TOB][i], color="blue", label=L"$\Delta$TOB1 = " * dTOB * " °C")
        fill_between(ctd[date][:t][i], ctd[date][:temp_TOB][i].+0.1, ctd[date][:temp_TOB][i].-0.1, color="blue", alpha=0.2)
        plot(ctd[date][:t][i],ctd[date][:temp_PT][i],color="red", label=L"$\Delta$T = " * dT * " °C")
        fill_between(ctd[date][:t][i], ctd[date][:temp_PT][i].+0.01, ctd[date][:temp_PT][i].-0.01, color="red", alpha=0.2)
        axhline(y=0,linewidth=1,color="black")
        ax = gca()
        ax.xaxis.set_major_formatter(majorformatter)
        if date == "2306"
            ax.xaxis.set_major_locator(majorlocator)
        end
        title("CTD-" * ctd["number"])
        if k == 1
            ylabel("Temperature [°C]")
        end
        xlabel("Hour of the day")
        grid(true)
        legend(handlelength=1.0)
    end
    gcf()
end

"""
Plot derived quantities (dphi/dz, Q, v, S, f, n') for all experiments.
Filled symbols correspond to experiments selected for further calculations and which are shown in Fig. 3 of the paper.
"""
function plot_data_unselected(date, ctd309, ctd265, ctd145, mid1, mid2, mid3, pick, error_p)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14
    fig, axs = subplots(6, 1, sharex=true, figsize=(16,9))
    fig.subplots_adjust(hspace=0.1) # Remove horizontal space between axes
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")
    pl1, pl2, pl3 = (nothing for i=1:100) # avoids introducing globals within the if statements
    dl309_265 = mid2[date][:dz]

    window_am15 = 60 # half of the window over which is averaged in the boxcar function
    window_am13 = 900

    # Hydraulic gradient
    idx = 1:minimum([length(ctd309[date][:t]), length(ctd265[date][:t])])
    X = map(Float64,convert(Array,idx))
    dphi = ((Particles.(n_partcl, Normal.(ctd265[date][:press][idx], error_p)) .-                              # hydraulic gradient in mH20/m
             Particles.(n_partcl, Normal.(ctd309[date][:press][idx], error_p)) ) ./  dl309_265 .- rhow*g) ./
            (rhow*g)
    if any(date .== ["2008", "2108"])
        window = window_am13 # smaller discharge requires different averaging window
    else
        window = window_am15
    end
    dphi_smooth = boxcar(pmean.(dphi), window)
    std_smooth = boxcar(pstd.(dphi), window)
    axs[1].plot(ctd309[date][:t][idx], dphi_smooth, "k", linewidth=0.5, label="CTD-309/CTD-265")
    #axs[1].fill_between(ctd309[date][:t][idx], dphi_smooth.+std_smooth, dphi_smooth.-std_smooth, color="grey", alpha=0.5)
    if haskey(ctd145,date)
        dl145_309 = mid3[date][:dz]
        idx = 1:minimum([length(ctd309[date][:t]), length(ctd145[date][:t])])
        X = map(Float64,convert(Array,idx))                                                                   # hydraulic gradient in mH20/m
        dphi = ((Particles.(n_partcl, Normal.(ctd309[date][:press][idx], error_p)) .-
                 Particles.(n_partcl, Normal.(ctd145[date][:press][idx], error_p)))   ./dl145_309 .- rhow*g) ./
               (rhow*g)
        dphi_smooth = boxcar(pmean.(dphi), window)
        std_smooth = boxcar(pstd.(dphi), window)
        axs[1].plot(ctd309[date][:t][idx], dphi_smooth, "-.k", linewidth=0.5, label="CTD-145/CTD-309")
        #axs[1].fill_between(ctd309[date][:t][idx], dphi_smooth.+std_smooth,dphi_smooth.-std_smooth, color="grey", alpha=0.5)
    end
    if haskey(ctd145, date)
        axs[1].legend(ncol=2)
    else
        axs[1].legend()
    end
    axs[1].set_ylabel(L"d\phi/dz\,\mathrm{[m\,H_2O/m]}")
    axs[1].xaxis.set_tick_params(bottom=false)
    axs[1].set_title(date[1:2] * "-Aug")

    # Discharge
    pl1 = axs[2].errorbar(mid1[date][:t_inj],pmean.(ctd309[date][:Q])*1e3,yerr=pstd.(ctd309[date][:Q])*1e3,fmt="o",color="darkblue",markerfacecolor="none",label="CTD-309, unreliable")
    pl2 = axs[2].errorbar(mid1[date][:t_inj],pmean.(ctd265[date][:Q])*1e3,yerr=pstd.(ctd265[date][:Q])*1e3,fmt="o",color="orange",markerfacecolor="none",label="CTD-265, unreliable")
    if haskey(pick, date)
        axs[2].errorbar(mid1[date][:t_inj][pick[date]],pmean.(ctd309[date][:Q][pick[date]])*1e3,yerr=pstd.(ctd309[date][:Q][pick[date]])*1e3,fmt="o",color="darkblue",label="CTD-309, reliable")
        axs[2].errorbar(mid1[date][:t_inj][pick[date]],pmean.(ctd265[date][:Q][pick[date]])*1e3,yerr=pstd.(ctd265[date][:Q][pick[date]])*1e3,fmt="o",color="orange",label="CTD-265, reliable")
    end
    if haskey(ctd145, date)
        pl3 = axs[2].errorbar(mid1[date][:t_inj],pmean.(ctd145[date][:Q])*1e3,yerr=pstd.(ctd145[date][:Q])*1e3,fmt="o",color="darkgreen",markerfacecolor="none",label="CTD-145, unreliable")
        if haskey(pick, date)
            axs[2].errorbar(mid1[date][:t_inj][pick[date]],pmean.(ctd145[date][:Q][pick[date]])*1e3,yerr=pstd.(ctd145[date][:Q][pick[date]])*1e3,fmt="o",color="darkgreen",label="CTD-145, reliable")
        end
    end
    axs[2].set_ylabel(L"$Q\,\mathrm{[l/s]}$")
    axs[2].xaxis.set_tick_params(bottom=false)
    if date == "2008"
        axs[2].legend((pl1,pl2,pl3),("CTD-309","CTD-265","CTD-145"),
                #loc="upper left"
                )
    elseif date == "2108"
        axs[2].legend((pl1,pl2,pl3),("CTD-309","CTD-265","CTD-145"),
                loc="upper center",
                ncol=3)
    else
        axs[2].legend((pl1,pl2),("CTD-309","CTD-265"),
                #loc="upper left"
                )
    end

    # Mean flow speed
    if haskey(ctd145,date)
        pl2 = axs[3].errorbar(mid2[date][:t_inj], pmean.(mid2[date][:v]), yerr=pstd.(mid2[date][:v]), fmt="ks", markerfacecolor="none")
        pl3 = axs[3].errorbar(mid3[date][:t_inj], pmean.(mid3[date][:v]), yerr=pstd.(mid3[date][:v]), fmt="kd", markerfacecolor="none")
        if haskey(pick, date)
            axs[3].errorbar(mid2[date][:t_inj][pick[date]], pmean.(mid2[date][:v][pick[date]]), yerr=pstd.(mid2[date][:v][pick[date]]), fmt="ks")
            axs[3].errorbar(mid3[date][:t_inj][pick[date]], pmean.(mid3[date][:v][pick[date]]), yerr=pstd.(mid3[date][:v][pick[date]]), fmt="kd")
        end
        if date == "2108"
            axs[3].legend((pl2,pl3),("surface/CTD-309","CTD-309/CTD-265"),
                loc="upper right",
                ncol=2
                )
        else
            axs[3].legend((pl2,pl3),("CTD-145/CTD-309","CTD-309/CTD-265"))
        end
    else
        pl1 = axs[3].errorbar(mid1[date][:t_inj],pmean.(mid1[date][:v]),yerr=pstd.(mid1[date][:v]),fmt="ko",markerfacecolor="none",label="CTD-309, unreliable")
        pl2 = axs[3].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:v]),yerr=pstd.(mid2[date][:v]),fmt="ks",markerfacecolor="none",label="CTD-265, unreliable")
        if haskey(pick,date)
            axs[3].errorbar(mid1[date][:t_inj][pick[date]],pmean.(mid1[date][:v][pick[date]]),yerr=pstd.(mid1[date][:v][pick[date]]),fmt="ko",label="CTD-309, reliable")
            axs[3].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:v][pick[date]]),yerr=pstd.(mid2[date][:v][pick[date]]),fmt="ks",label="CTD-265, reliable")
        end
        axs[3].legend((pl1,pl2),("surface/CTD-309","CTD-309/CTD-265"))
    end
    axs[3].set_ylabel(L"$v\,\mathrm{[m/s]}$")
    axs[3].xaxis.set_tick_params(bottom=false)

    # Cross-sectional area
    if haskey(ctd145,date)
        axs[4].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:S]),yerr=pstd.(mid2[date][:S]),fmt="ks",markerfacecolor="none",label="CTD-265, unreliable")
        axs[4].errorbar(mid3[date][:t_inj],pmean.(mid3[date][:S]),yerr=pstd.(mid3[date][:S]),fmt="kd",markerfacecolor="none",label="CTD-145, unreliable")
        if haskey(pick, date)
            axs[4].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:S][pick[date]]),yerr=pstd.(mid2[date][:S][pick[date]]),fmt="ks",label="CTD-265, reliable")
            axs[4].errorbar(mid3[date][:t_inj][pick[date]],pmean.(mid3[date][:S][pick[date]]),yerr=pstd.(mid3[date][:S][pick[date]]),fmt="kd",label="CTD-145, reliable")
        end
    else
        axs[4].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:S]),yerr=pstd.(mid2[date][:S]),fmt="ks",markerfacecolor="none",label="CTD-309/CTD-265")
        if haskey(pick,date)
            axs[4].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:S][pick[date]]),yerr=pstd.(mid2[date][:S][pick[date]]),fmt="ks",label="")
        end
    end
    axs[4].set_ylabel(L"$S\,\mathrm{[m^2]}$")
    axs[4].xaxis.set_tick_params(bottom=false)

    # friction factor
    if haskey(ctd145,date)
        axs[5].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:f]),yerr=pstd.(mid2[date][:f]),fmt="ks",markerfacecolor="none",label="CTD-265, unreliable")
        axs[5].errorbar(mid3[date][:t_inj],pmean.(mid3[date][:f]),yerr=pstd.(mid3[date][:f]),fmt="kd",markerfacecolor="none",label="CTD-145, unreliable")
        if haskey(pick, date)
            axs[5].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:f][pick[date]]),yerr=pstd.(mid2[date][:f][pick[date]]),fmt="ks",label="CTD-265, reliable")
            axs[5].errorbar(mid3[date][:t_inj][pick[date]],pmean.(mid3[date][:f][pick[date]]),yerr=pstd.(mid3[date][:f][pick[date]]),fmt="kd",label="CTD-145, reliable")
        end
    else
        axs[5].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:f]),yerr=pstd.(mid2[date][:f]),fmt="ks",markerfacecolor="none",label="CTD-309/CTD-265")
        if haskey(pick,date)
            axs[5].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:f][pick[date]]),yerr=pstd.(mid2[date][:f][pick[date]]),fmt="ks",label="")
        end
        #axs[5].set_yscale("log")
    end
    axs[5].set_ylabel(L"$f$")
    axs[5].xaxis.set_tick_params(bottom=false)

    # manning roughness
    if haskey(ctd145,date)
        axs[6].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:n_manning]),yerr=pstd.(mid2[date][:n_manning]),fmt="ks",markerfacecolor="none",label="CTD-265, unreliable")
        axs[6].errorbar(mid3[date][:t_inj],pmean.(mid3[date][:n_manning]),yerr=pstd.(mid3[date][:n_manning]),fmt="kd",markerfacecolor="none",label="CTD-145, unreliable")
        if haskey(pick, date)
            axs[6].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:n_manning][pick[date]]),yerr=pstd.(mid2[date][:n_manning][pick[date]]),fmt="ks",label="CTD-265, reliable")
            axs[6].errorbar(mid3[date][:t_inj][pick[date]],pmean.(mid3[date][:n_manning][pick[date]]),yerr=pstd.(mid3[date][:n_manning][pick[date]]),fmt="kd",label="CTD-145, reliable")
        end
    else
        axs[6].errorbar(mid2[date][:t_inj],pmean.(mid2[date][:n_manning]),yerr=pstd.(mid2[date][:n_manning]),fmt="ks",markerfacecolor="none",label="CTD-309/CTD-265")
        if haskey(pick,date)
            axs[6].errorbar(mid2[date][:t_inj][pick[date]],pmean.(mid2[date][:n_manning][pick[date]]),yerr=pstd.(mid2[date][:n_manning][pick[date]]),fmt="ks",label="")
        end
        #axs[6].set_yscale("log")
    end
    axs[6].set_ylabel(L"$n'$")
    axs[6].xaxis.set_major_formatter(majorformatter)
    axs[6].set_xlabel("Time")

    gcf()
end

"""
Plot the tracer signals at three sensors (only works 21 August)
"""
function plot_example_traces(ctds, indices, date, tr)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 8 # general
    formatter = matplotlib.dates.DateFormatter("%M")
    msec_per_year = 3600*24*365*1e3

    figure(figsize=(6,3), dpi=pixel_resolution)
    ax=PyPlot.axes()
    for (ctd, idx, col) in  zip(ctds,
                    [indices[date][5][tr]:indices[date][6][tr],     # indices of ctd-145
                    indices[date][1][tr]:indices[date][2][tr]-1700, # indices of ctd-309
                    indices[date][3][tr]:indices[date][4][tr]-1700], # indices of ctd-265
                    ["darkblue", "orange", "darkgreen"])
        dt = convert.(DateTime,ctd[date][:t][idx].-ctd[date][:t][idx[1]].+convert(DateTime,Dates.Millisecond(msec_per_year*1))) # PyPlot doesn't want to plot year zero
        plot(dt, ctd[date][:cond][idx], label= "CTD-" * ctd["number"], color=col)
    end
    ax.xaxis.set_major_formatter(formatter)
    xlabel("Minutes after injection")
    ylabel("Conductivity [μS/cm]")
    legend()
    gcf()
end

"""
Plot the raw CTD measurements: conductivity, temperature (already corrected using zero degree calibration) and pressure
"""
function plot_CTD(date, ctd309, ctd265, ctd145, e_T, range)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12 # general
    fig, axs = subplots(3, 1, sharex=true, figsize=(12,7))
    fig.subplots_adjust(hspace=0) # Remove horizontal space between axes
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")
    range_dpress = floor.(Int, range./120)
    if date=="2108"
        range_dpress = 1:length(ctd309[date][:dpress])
    end

    # conductivity
    axs[1].plot(ctd309[date][:t][range],ctd309[date][:cond][range], label="CTD-309",color="darkblue")
    axs[1].plot(ctd265[date][:t][range],ctd265[date][:cond][range], label="CTD-265",color="orange")
    if haskey(ctd145,date)
        axs[1].plot(ctd145[date][:t][range],ctd145[date][:cond][range], label="CTD-145",color="darkgreen")
    end
    axs[1].legend(loc="upper left")
    axs[1].set_ylabel("Conductivity [μS/cm]")
    axs[1].set_title(date[1:2] * "-Aug")
    if date=="0808"
        axs[1].set_ylim(0,400)
    end

    # pressure
    axs[2].plot(ctd309[date][:t][range],ctd309[date][:press][range] ./ (rhow*g), label="CTD-309",color="darkblue")
    axs[2].plot(ctd265[date][:t][range],ctd265[date][:press][range] ./ (rhow*g), label="CTD-265",color="orange")
    if haskey(ctd145,date)
        axs[2].plot(ctd145[date][:t][range],ctd145[date][:press][range] ./ (rhow*g), label="CTD-145",color="darkgreen")
    end
    axs[2].legend(loc="upper left")
    axs[2].set_ylabel(L"\mathrm{Pressure\,[m\,H_2O]}")

    # temperature
    axs[3].plot(ctd309[date][:t][range], ctd309[date][:temp_PT][range], label="CTD-309",color="darkblue")
    axs[3].fill_between(ctd309[date][:t][range], ctd309[date][:temp_PT][range].-e_T, ctd309[date][:temp_PT][range].+e_T, color="darkblue", alpha=0.3, zorder=1)
    #axs[3].fill_between(ctd309[date][:t_dpress][range_dpress], pmean.(ctd309[date][:dpress][range_dpress]*ct)-pstd.(ctd309[date][:dpress][range_dpress]*ct), pmean.(ctd309[date][:dpress][range_dpress]*ct)+pstd.(ctd309[date][:dpress][range_dpress]*ct), color="blue", zorder=2, label="CTD-309, T_melt range")
    axs[3].plot(ctd265[date][:t][range], ctd265[date][:temp_PT][range], label="CTD-265",color="orange")
    axs[3].fill_between(ctd265[date][:t][range], ctd265[date][:temp_PT][range].-e_T, ctd265[date][:temp_PT][range].+e_T, color="orange", alpha=0.3, zorder=1)
    #axs[3].fill_between(ctd265[date][:t_dpress][range_dpress], pmean.(ctd265[date][:dpress][range_dpress]*ct)-pstd.(ctd265[date][:dpress][range_dpress]*ct), pmean.(ctd265[date][:dpress][range_dpress]*ct)+pstd.(ctd265[date][:dpress][range_dpress]*ct), color="red", zorder=2, label="CTD-309, T_melt range")
    if haskey(ctd145,date)
        axs[3].plot(ctd145[date][:t][range],ctd145[date][:temp_PT][range], label="CTD-145",color="darkgreen")
        axs[3].fill_between(ctd145[date][:t][range], ctd145[date][:temp_PT][range].-e_T, ctd145[date][:temp_PT][range].+e_T, color="darkgreen", alpha=0.3, zorder=1)
    end
    axs[3].legend(loc="upper left")
    axs[3].xaxis.set_major_formatter(majorformatter)
    axs[3].set_ylabel("Temperature [°C]")
    axs[3].set_xlabel("Hour of the day")
    gcf()
end

"""
Plot temperature together with pressure melting temperature
"""
function plot_temp(date, ctd309, ctd265, range, error_T)
    figure(figsize=(12,6), dpi=pixel_resolution)
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")
    ax1=subplot(211)
    range_dpress = floor.(Int, range./120)
    if date=="2108"
        range_dpress = 1:length(ctd309[date][:dpress])
    end
    plot(ctd309[date][:t][range], ctd309[date][:temp_PT][range], label="CTD-309", color="darkblue")
    fill_between(ctd309[date][:t][range], ctd309[date][:temp_PT][range].-error_T, ctd309[date][:temp_PT][range].+error_T, color="darkblue", alpha=0.3, zorder=1)
    fill_between(ctd309[date][:t_dpress][range_dpress], ctd309[date][:dpress][range_dpress]*ct_pure, ctd309[date][:dpress][range_dpress]*ct_air, color="blue", zorder=2, label="CTD-309, T_melt range") # error of p not so important
    legend()
    ax1.xaxis.set_major_formatter(majorformatter)
    ylabel("Temperature [°C]")
    title(date[1:2] * "-Aug")

    ax2=subplot(212, sharex=ax1)
    plot(ctd265[date][:t][range], ctd265[date][:temp_PT][range], label="CTD-265", color="orange")
    fill_between(ctd265[date][:t][range], ctd265[date][:temp_PT][range].-error_T, ctd265[date][:temp_PT][range].+error_T, color="orange", alpha=0.3, zorder=1)
    fill_between(ctd265[date][:t_dpress][range_dpress], ctd265[date][:dpress][range_dpress]*ct_pure, ctd265[date][:dpress][range_dpress]*ct_air, color="red", zorder=2, label="CTD-265, T_melt range")
    legend()
    ax2.xaxis.set_major_formatter(majorformatter)
    ylabel("Temperature [°C]")
    xlabel("Hour of the day")
    gcf()
end

"""
Plot only closure rates of 09-Aug and 21-Aug
"""
function plot_closure(mid_309_265, model_runs)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14 # general
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")

    figure(figsize=(8,5), dpi=pixel_resolution)
    ax = subplot(1, 1, 1)
    errorbar(mid_309_265["0908"][:t_inj] , pmean.(model_runs["0908"][:closure]),yerr=pstd.(model_runs["0908"][:closure]), fmt="k_", mew=2, ms=8, label="9-Aug/AM15, positive")
    errorbar(mid_309_265["2108"][:t_inj] .- Day(12) ,abs.(pmean.(model_runs["2108"][:closure])),yerr=pstd.(model_runs["2108"][:closure]), fmt="g_", mew=2, ms=8, label="21-Aug/AM13, negative")
    yscale("log")
    xlabel("Time of the corresponding day")
    ylabel(L"$|v_c|\,[\mathrm{m^2/s}]$")
    legend()
    ax.xaxis.set_major_formatter(majorformatter)
    gcf()
end


"""
Plot opening rates of ct-gradient and free-gradient model, as well as the frictional component which is the same for both
"""
function plot_opening(mid_309_265, model_runs)

    output = Dict("0908" => Dict(),
                  "2108" => Dict())
    for date in ["0908", "2108"]
        # the results from the MCMC runs need to be averaged first to obtain the most likely time series
        output[date][:dSdt_MCMC] = Particles.(n_partcl, Normal.(mean(model_runs[date][:dSdt_MCMC]), std(model_runs[date][:dSdt_MCMC])))
        output[date][:dSdt_ct] = model_runs[date][:dSdt_ct]
        output[date][:dSdt_frictional] = model_runs[date][:dSdt_ct] .- model_runs[date][:dSdt_sensible_ct]
    end

    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14 # general
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")
    # majorlocator = matplotlib.dates.HourLocator(byhour=(11, 13, 15, 17))
    figure(figsize=(16,5), dpi=pixel_resolution)

    for (nd, (date, tit, lab)) in enumerate(zip(["0908", "2108"], ["9-Aug/AM15", "21-Aug/AM13"], [L"\bf{a}", L"\bf{b}"]))
        ax = subplot(1, 2, nd)

        cols = ["black", "grey", "royalblue"]
        syms = ["_", "_", "o"]
        symsize = [9, 9, 5]
        labs = ["ct-gradient model, total",
                 "Free-gradient model, total",
                 "Frictional heat component"]
        mods = [:dSdt_ct, :dSdt_MCMC, :dSdt_frictional ]

        for p in 1:length(cols)
            errorbar(mid_309_265[date][:t_inj], pmean.(output[date][mods[p]]), yerr=pstd.(output[date][mods[p]]), color=cols[p], fmt=syms[p], ms=symsize[p], label=labs[p])
        end

        ax.xaxis.set_major_formatter(majorformatter)


        if date == "0908"
            ylabel(L"Opening rate [$\mathrm{m^2/s}$]")
        end
        if date == "2108"
            legend()
        end

        xlabel("Time")
        title(tit)
        text(-0.1, 1.05, lab, transform=ax.transAxes, ha="left", va="bottom", fontsize=19)
    end
    gcf()
end

"""
Plot opening rates of ct-model and free-gradient model vs. closure rates to compare their magnitudes
"""
function plot_opening_closure(mid_309_265, model_runs)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14 # general
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")
    # majorlocator = matplotlib.dates.HourLocator(byhour=(11, 13, 15, 17))
    figure(figsize=(16,5), dpi=pixel_resolution)

    a1, a2, a3, a4 = (nothing for i=1:100)
    for (nd, (date, lab)) in enumerate(zip(["0908", "2108"], [L"\bf{a}", L"\bf{b}"]))
        ax = subplot(1, 2, nd)
        #ax1=PyPlot.axes([0.1, 0.1, 0.4, 0.8])  # left, bottom, width, height
        a1 = errorbar(mid_309_265[date][:t_inj],pmean.(model_runs[date][:dSdt_ct]), yerr=pstd.(model_runs[date][:dSdt_ct]),fmt="k.", ms=8)
        a2 = errorbar(mid_309_265[date][:t_inj],mean(model_runs[date][:dSdt_MCMC]), yerr=std(model_runs[date][:dSdt_MCMC]),fmt="b.", ms=8)
        if date == "0908"
            a3 = errorbar(mid_309_265[date][:t_inj],pmean.(model_runs[date][:closure]),yerr=pstd.(model_runs[date][:closure]), fmt="g+", mew=2, ms=8)
        else
            a4 = errorbar(mid_309_265[date][:t_inj],abs.(pmean.(model_runs[date][:closure])),yerr=pstd.(model_runs[date][:closure]), fmt="gx", mew=2, ms=8)
        end
        ax.xaxis.set_major_formatter(majorformatter)
        # ax.xaxis.set_major_locator(majorlocator)
        yscale("log")
        text(0.0, 1.0, lab, transform=ax.transAxes, ha="left", va="bottom", fontsize=19)
        xlabel("Time")
        if date == "0908"
            ylabel(L"Closure/Opening rate [$\mathrm{m^2/s}$]")
        else
            legend((a1, a2, a3, a4),
                (L"opening $c_{t}$-model",
                 "opening free-gradient model",
                "creep closure positive",
                "creep closure netagive"),
                ncol = 1
                )
        end
        title("AM15, " * date[1:2] * "-Aug")
    end
    gcf()
end

"""
Plot for selected data points: hydraulic gradient, discharge, velocity, cross-sectional area, friction factor and manning roughness.
Fig. 3 in the paper.
"""
function multiplot(mid_309_265, pick, ctd309, ctd265, e_p, idx_plot, idx_gaps)
    dates = ["0808", "0908", "1008", "1108", "1308", "2108"]
    props = [:dphi_dz, :Q, :v, :S, :Re, :f, :n_manning]
    panel_labs = [L"\bf{a}", L"\bf{b}", L"\bf{c}", L"\bf{d}", L"\bf{e}", L"\bf{f}", L"\bf{g}"]
    ylabels = [#L"$p_w\,\mathrm{(mH_2O)}$",
               L"$\partial\phi/\partial z$ $\mathrm{(mH_2O\,m^{-1})}$",
               L"$Q$ $\mathrm{(m^3\,s^{-1})}$", # or in l/s ???
               L"$v$ $\mathrm{(m\,s^{-1})}$",
               L"$S$ $\mathrm{(m^2)}$",
               L"$Re$",
               L"$f$",
               L"$n'$ $(\mathrm{s\,m^{-1/3}})$"]
    nprops = length(props)

    # font properties
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 20
    fs = 20 # font size for figure numbering

    # time axis format
    # for the Locator it is necessary to define this function,
    # otherwise the xticks from previous plots are removed
    # https://stackoverflow.com/questions/55010236/set-major-locator-removes-x-ticks-and-labels-from-previous-subplots

    function format_xaxis(ax, loc)
        fmt = matplotlib.dates.DateFormatter("%H:%M")
        loc = matplotlib.dates.HourLocator(;loc...)
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(loc)
    end

    # width ratios for left subfigure
    widths = []
    for date in dates[1:end-1]
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    w_space = 0.04
    grid_dict_left = Dict(:width_ratios => widths,
                          :wspace => w_space, # horizontal space between panels
                          :hspace => 0.)  # vertical space between panels

    # grid_dict for right subfigure
    grid_dict_right = Dict(:hspace => 0.)

    # style of the spine line at x-axis breakpoints
    spine_line = (0, (2, 4))

    # for drawing breakpoints
    break_point_length = 0.017
    dxs = break_point_length ./ w_ratios
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :linewidth => 1.0,
                  :clip_on   => false
                  )

    # draw subplots
    fig = plt.figure(figsize=(20,15), dpi=pixel_resolution)
    subfigs = fig.subfigures(1, 2, width_ratios=[sum(widths)+mean(widths)*4*w_space, length(idx_plot["2108"])], wspace = 0.0)
    axesleft  = subfigs[1].subplots(length(props), 5, sharey="row", sharex="col", gridspec_kw=grid_dict_left)
    axesright = subfigs[2].subplots(length(props), 1, sharex="col", gridspec_kw=grid_dict_right)

    for (nd, date) in enumerate(dates)
        i = idx_plot[date]

        # time axes injections
        t_inj = mid_309_265[date][:t_inj][pick[date]]

        # pressure
        # time windows over which to average pressures, in seconds
        if date == "2108"
            window = 900
        else
            window = 60
        end
        # smooth the pressure difference
        presstop = boxcar(ctd309[date][:press][i], window)
        pressbot = boxcar(ctd265[date][:press][i], window)

        # hydraulic gradient
        dz = mid_309_265[date][:dz]
        dpress = Particles.(n_partcl, Normal.(pressbot, e_p*300)) .-
                 Particles.(n_partcl, Normal.(presstop, e_p*100))
        #dp_smooth = Particles.(n_partcl, Normal.(pmean.(dpress), pstd.(dpress)))
        dphi_dz = dpress ./ dz .- rhow*g # compute hydraulic gradient, Pa/m

         # conversion to mH2O
        presstop = presstop ./ (rhow * g)
        pressbot = pressbot ./ (rhow * g)
        dphi_dz = dphi_dz ./ (rhow*g)

        # remove data points where we moved CTDs up and down or blocked the water inflow
        if haskey(idx_gaps, date)
            presstop[idx_gaps[date]] .= NaN
            pressbot[idx_gaps[date]] .= NaN
            dphi_dz[idx_gaps[date]] .= NaN
        end

        for (row, (prop, ylab)) in enumerate(zip(props, ylabels))
            if date == "2108"
                ax = axesright[row]
                nd = 1
            else
                ax = axesleft[nprops*(nd-1) + row]
            end
            if row == 1 # first row is for hydraulic gradient, manually
                ax.plot(ctd309[date][:t][i], pmean.(dphi_dz), "k", linewidth=0.5)
                ax.fill_between(ctd309[date][:t][i], pmean.(dphi_dz) .+ pstd.(dphi_dz), pmean.(dphi_dz) .- pstd.(dphi_dz), color="grey")
            else
                data = mid_309_265[date][prop][pick[date]]
                ax.errorbar(t_inj, pmean.(data), yerr=pstd.(data), fmt="kx", markersize=7)
            end

            # logarithmic y-scale for some parameters
            if any(prop .== [:f, :n_manning, :Re])
                ax.set_yscale("log")
            end

            # title
            if row == 1
                ax.set_title(date[1:2] * "-Aug")
            end

            # y-label and panel label
            if date == "0808"
                ax.set_ylabel(ylab, labelpad=25., wrap=true) #, rotation="horizontal", va="center")
                text(0.1, 0.6, panel_labs[row], fontsize=fs, transform=ax.transAxes, ha="left", va="top")
            end

            # make intermediate spines dashed lines and remove axis ticks
            if date !== "2108"
                dx = dxs[nd]
                dy = break_point_length * 3 # needs the factor because plots longer than high, factor chosen randomly
                if nd !== length(dates)-1 # adjust right side of panel
                    ax.spines["right"].set_linestyle(spine_line)
                    ax.tick_params(right=false)
                end
                if nd !== 1 # adjust left side of panel
                    ax.spines["left"].set_linestyle(spine_line)
                    ax.tick_params(which="both", left=false) # both -> remove also minor ticks in logscale plots, default is major
                end
            end
            if row !== nprops
                ax.tick_params(bottom=false)
            end

            # draw breakpoints
            kwargs[:transform] = ax.transAxes
            if row == nprops && date !== "2108"
                if nd !== 1 # draw in lower left corner
                    ax.plot((-dx,+dx), (-dy,+dy); kwargs...)
                end
                if nd !== length(dates)-1 # draw in lower right corner
                    ax.plot((1-dx,1+dx),(-dy,+dy); kwargs...)
                end
            end
            if date == "1308"
                format_xaxis(ax, Dict(:byhour=>14))
            else
                format_xaxis(ax, Dict(:interval=>2))
            end
            ax.tick_params("x", pad=10., labelrotation=30.)

            if prop == :Re && date != "2108"
                ax.margins(y=0.5) # so that ylabels don't overlap each other
            elseif prop == :n_manning && date != "2108"
                ax.margins(y=0.65)
            else
                ax.margins(y=0.2)
            end
        end
    end

    subfigs_pos = Dict(:left => 0.18, :right => 0.97, :bottom => 0.12, :top => 0.92)
    subfigs[1].subplots_adjust(;subfigs_pos...)
    subfigs[2].subplots_adjust(;subfigs_pos...)
    subfigs[1].suptitle(L"\bf{AM15}", y=0.98)
    subfigs[2].suptitle(L"\bf{AM13}", y=0.98)
    subfigs[1].supxlabel("Time", y=0.02)

    gcf()

end

"""
Produce Fig. S1 (supplements) showing water pressure and temperature
"""
function plot_pw_Tw_supp(mid_309_265, pick, ctd309, ctd265, e_p, e_T, idx_plot, idx_gaps)
    dates = ["0808", "0908", "1008", "1108", "1308", "2108"]
    props = [:press, :temp_PT]
    panel_labs = [L"\bf{a}", L"\bf{b}"]
    ylabels = [L"$p_w\,\mathrm{(mH_2O)}$",
               L"$T_w\,(\degree\,\mathrm{C})$"]
    errors = [[100*e_p, 300*e_p], [e_T, e_T]]
    nprops = length(props)

    # font properties
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14
    fs = 14 # font size for figure numbering

    # legend handles
    L2D = PyPlot.matplotlib.lines.Line2D
    Patch = PyPlot.matplotlib.patches.Patch
    custom_legend = [L2D([0], [0], color="darkblue", lw=2),
                     L2D([0], [0], color="orange", lw=2),
                     Patch(facecolor="darkblue", alpha = 0.3),
                     Patch(facecolor="orange", alpha = 0.3)]

    # time axis format
    # for the Locator it is necessary to define this function,
    # otherwise the xticks from previous plots are removed
    # https://stackoverflow.com/questions/55010236/set-major-locator-removes-x-ticks-and-labels-from-previous-subplots

    function format_xaxis(ax, loc)
        fmt = matplotlib.dates.DateFormatter("%H:%M")
        loc = matplotlib.dates.HourLocator(;loc...)
        ax.xaxis.set_major_formatter(fmt)
        ax.xaxis.set_major_locator(loc)
    end

    # width ratios for left subfigure
    widths = []
    for date in dates[1:end-1]
        append!(widths, length(idx_plot[date]))
    end
    w_ratios = widths ./ mean(widths)
    w_space = 0.04
    grid_dict_left = Dict(:width_ratios => widths,
                          :wspace => w_space, # horizontal space between panels
                          :hspace => 0.)  # vertical space between panels

    # grid_dict for right subfigure
    grid_dict_right = Dict(:hspace => 0.)

    # style of the spine line at x-axis breakpoints
    spine_line = (0, (2, 4))

    # for drawing breakpoints
    break_point_length = 0.01
    dxs = break_point_length ./ w_ratios
    kwargs = Dict(:transform => nothing,
                  :color     => "k",
                  :linewidth => 1.0,
                  :clip_on   => false
                  )

    # draw subplots
    fig = plt.figure(figsize=(14,6), dpi=pixel_resolution)
    subfigs = fig.subfigures(1, 2, width_ratios=[sum(widths)+mean(widths)*4*w_space, length(idx_plot["2108"])], wspace = 0.0)
    axesleft  = subfigs[1].subplots(length(props), 5, sharey="row", sharex="col", gridspec_kw=grid_dict_left)
    axesright = subfigs[2].subplots(length(props), 1, sharex="col", gridspec_kw=grid_dict_right)

    # share y-axis between AM15 and AM13
    for a in 1:length(axesright)
        axesright[a].sharey(axesleft[(a-1)*5+1])
    end
    for (nd, date) in enumerate(dates)
        i = idx_plot[date]

        # time axes injections
        t_inj = mid_309_265[date][:t_inj][pick[date]]

        # pressure
        # time windows over which to average pressures, in seconds
        if date == "2108"
            window = 900
        else
            window = 60
        end

        for (row, (prop, ylab, err)) in enumerate(zip(props, ylabels, errors))

            top = boxcar(ctd309[date][prop][i], window)
            bot = boxcar(ctd265[date][prop][i], window)
            # remove data points where we moved CTDs up and down or blocked the water inflow
            if haskey(idx_gaps, date)
                top[idx_gaps[date]] .= NaN
                bot[idx_gaps[date]] .= NaN
            end

            top = Particles.(n_partcl, Normal.(top, err[1]))
            bot = Particles.(n_partcl, Normal.(bot, err[2]))

            if prop == :press # pressure
                # conversion to mH2O
                top = top ./ (rhow * g)
                bot = bot ./ (rhow * g)
            end

            if date == "2108"
                ax = axesright[row]
                nd = 1
            else
                ax = axesleft[nprops*(nd-1) + row]
            end

            ax.plot(ctd309[date][:t][i], pmean.(top), color="darkblue", linewidth=2, label="upper sensor")
            ax.plot(ctd309[date][:t][i], pmean.(bot), color="orange", linewidth=2, label="lower sensor")
            ax.fill_between(ctd309[date][:t][i], pmean.(top) .+ pstd.(top), pmean.(top) .- pstd.(top), color="darkblue", alpha=0.3, label="uncertainty")
            ax.fill_between(ctd309[date][:t][i], pmean.(bot) .+ pstd.(bot), pmean.(bot) .- pstd.(bot), color="orange", alpha=0.3, label="uncertainty")
            legend(custom_legend,
                   ("Upper sensor", "Lower sensor", "Uncertainty", "Uncertainty"),
                   ncol = 2,
                   loc = "upper center",
                   handlelength = 0.8,
                   columnspacing = 1.0,
                   frameon = false)
            # title
            if row == 1
                ax.set_title(date[1:2] * "-Aug")
            end

            # y-label and panel label
            if date == "0808"
                ax.set_ylabel(ylab, labelpad=25., wrap=true) #, rotation="horizontal", va="center")
                text(0.1, 0.9, panel_labs[row], fontsize=fs, transform=ax.transAxes, ha="left", va="top")
            end

            # make intermediate spines dashed lines and remove axis ticks
            if date !== "2108"
                dx = dxs[nd]
                dy = break_point_length * 2 # needs the factor because plots longer than high, factor chosen randomly
                if nd !== length(dates)-1 # adjust right side of panel
                    ax.spines["right"].set_linestyle(spine_line)
                    ax.tick_params(right=false)
                end
                if nd !== 1 # adjust left side of panel
                    ax.spines["left"].set_linestyle(spine_line)
                    ax.tick_params(which="both", left=false) # both -> remove also minor ticks in logscale plots, default is major
                end
            end
            if row !== nprops
                ax.tick_params(bottom=false)
            end

            # draw breakpoints
            kwargs[:transform] = ax.transAxes
            if row == nprops && date !== "2108"
                if nd !== 1 # draw in lower left corner
                    ax.plot((-dx,+dx), (-dy,+dy); kwargs...)
                end
                if nd !== length(dates)-1 # draw in lower right corner
                    ax.plot((1-dx,1+dx),(-dy,+dy); kwargs...)
                end
            end
            # adjust format of time axis
            if date == "1308"
                format_xaxis(ax, Dict(:byhour=>14))
            else
                format_xaxis(ax, Dict(:interval=>2))
            end
            ax.tick_params("x", pad=10., labelrotation=30.)
            ax.margins(y=0.2) # so that ylabels don't overlap each other
        end
    end

    subfigs_pos = Dict(:left => 0.18, :right => 0.97, :bottom => 0.2, :top => 0.88)
    subfigs[1].subplots_adjust(;subfigs_pos...)
    subfigs[2].subplots_adjust(;subfigs_pos...)
    subfigs[1].suptitle(L"\bf{AM15}", y=0.98)
    subfigs[2].suptitle(L"\bf{AM13}", y=0.98)
    subfigs[1].supxlabel("Time", y=0.02)

gcf()

end

"""
Plot evolution of cross-sectional area according to measurements, ct-model and free-gradient model;
Plot corresponding distributions of temperature gradient
Fig. 4 in the paper.
"""
function plot_model_outputs(mid_309_265, model_runs)
    # some numbers for drawing the plots
    global left = 0.08
    global bottomlow = 0.1
    global bottomup = 0.45
    global width = 0.4
    global heightlow = 0.25
    global heightup = 0.5

    # workaround to make legend for 2D lines
    L2D = PyPlot.matplotlib.lines.Line2D
    Patch = PyPlot.matplotlib.patches.Patch
    custom_legend = [L2D([0], [0], color="black", lw=0.4),
                     L2D([0], [0], color="black", lw=1),
                     Patch(facecolor="grey", alpha = 0.5),
                     Patch(facecolor="grey", alpha = 1.0),
                     Patch(facecolor=(0, 0, 1, 0.1), edgecolor=(0, 0, 1, 1))]

    # time axis format
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")

    # font size
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 15
    fs = 18 # font size for figure numbering

    figure(figsize=(17,8), dpi=pixel_resolution)

    ax, pl1, pl2 = (nothing for i=1:100)
    for (date, tit) in zip(["0908", "2108"], ["9-Aug/AM15", "21-Aug/AM13"])

        if date == "2108"
            global left = 0.95-width
        end

        # time
        t = mid_309_265[date][:t_inj]

        # temperature gradients
        dTdz_MCMC = model_runs[date][:dTdz_MCMC]
        dTdz_ct = model_runs[date][:dTdz_ct]
        dTdz_measured = model_runs[date][:dTdz_measured]

        # size evolution
        S_MCMC = model_runs[date][:S_MCMC]
        S_ct = model_runs[date][:S_ct]

        # plot model vs. data
        ax = PyPlot.axes([left, bottomup, width, heightup])
        pl1 = errorbar(t, pmean.(mid_309_265[date][:S]), yerr=pstd.(mid_309_265[date][:S]), fmt="k+", zorder=3, markersize=10) # data
        [plot(t, S_MCMC[i], "grey", linewidth=0.3, zorder=1) for i in 1:5000:length(S_MCMC)] # model outputs
        plot(t, pmean.(S_ct)-pstd.(S_ct), "b", lw=1, zorder=2)
        plot(t, pmean.(S_ct)+pstd.(S_ct), "b", lw=1, zorder=2)
        pl2 = fill_between(t, pmean.(S_ct)-pstd.(S_ct), pmean.(S_ct)+pstd.(S_ct), color="blue", alpha=0.1, zorder=2, label=L"Range for $c_t$ and initial $S$")
        ax.xaxis.set_major_formatter(majorformatter)
        ylabel(L"$S\,\mathrm{(m^2)}$")
        xlabel("Time")
        title(tit)
        if date == "0908"
            legend((pl1, custom_legend[5], custom_legend[1]),
                   ("Measurements", L"$c_t$-gradient model", "Free-gradient model"),
                   framealpha = 0.8,
                   loc = "lower right"
                   #bbox_to_anchor=(1.2, 0.5, 0.6, .102), # position of the legend (x0, y0, width, height)
            )
        end
        if date == "0908"
            text(0.05, 0.93, L"\bf{a}", fontsize=fs, transform=ax.transAxes, ha="left", va="top")
        else
            text(0.05, 0.93, L"\bf{b}", fontsize=fs, transform=ax.transAxes, ha="left", va="top")
        end

        # plot distributions of dTw/dz
        ax = PyPlot.axes([left, bottomlow, width, heightlow])
        hist(dTdz_MCMC, bins=100, color="grey", alpha=0.5, zorder=2, density=true)
        hist(dTdz_MCMC, bins=100, range=(quantile(dTdz_MCMC, 0.025), quantile(dTdz_MCMC, 0.975)), color="grey", alpha=1.0, density=true)
        xlabel(L"$\partial T_w/\partial z\,(\mathrm{°C\,m^{-1}})$")
        ylabel("Frequency")
        tick_params("y", direction="in", color="white", left=false, labelleft=false) # bit of a workaround, for some reason it doesn't respond to the false

        ax.twinx() # Create another y-axis on top of the current one
        if date == "0908"
            nbin = 35
        else
            nbin = 25
        end
        count, bins, ignored = hist(Array(dTdz_measured), bins=nbin, alpha=0.0) #, density=true, histtype="step", color="black")
        plot(bins[2:end], count, color="black") #, linewidth=2, color='r')
        tick_params("y", right=false, labelright=false, left=false, labelleft=false)

        ax.twinx()
        hist(Array(dTdz_ct), bins=12, color="blue", alpha=0.1)
        hist(Array(dTdz_ct), bins=12, color="blue", histtype="step")

        if date == "2108"
            legend((custom_legend[2], custom_legend[5], custom_legend[3], custom_legend[4]),
                   ("Measurements", L"$c_t$-gradient model", "Free-gradient model", "95% CI"),
                   loc = "lower right",
                   framealpha=0.8
                   #bbox_to_anchor=(1.11, 0.7, 0.5, .102)) # position of the legend (x0, y0, width, height
                   )
        end
        xlim([-0.0025, 0.001])
        tick_params("y", right=false, labelright=false, left=false, labelleft=false)

        if date == "0908"
            text(0.05, 0.86, L"\bf{c}", fontsize=fs, transform=ax.transAxes, ha="left", va="top")
        else
            text(0.05, 0.86, L"\bf{d}", fontsize=fs, transform=ax.transAxes, ha="left", va="top")
        end
    end

    gcf()
end

"""
Plot evolution of cross-sectional area together with mean linear opening rate from ct-model. Similar to plot shown at SGM.
"""
function plot_dSdt_linear(mid, model_runs)

    # format time axis
    majorformatter = matplotlib.dates.DateFormatter("%H:%M")

    # font properties
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 15
    rcParams["font.serif"] = "Times New Roman"

    # workaround to make legend for 2D lines
    L2D = PyPlot.matplotlib.lines.Line2D
    custom_lines = [L2D([0], [0], color="k", linestyle="-")]

    # draw figure
    figure(figsize=(18,7), dpi=pixel_resolution)
    pl1, pl2, pl3 = (nothing for i=1:100)
    for (d, (date, midpt, lab)) in enumerate(zip(["0908", "2108"], [7, 9], [L"\bf{a}", L"\bf{b}"]))
        ax=subplot(1,2,d)

        ## S measurements
        pl1 = errorbar(mid[date][:t_inj],pmean.(mid[date][:S]),yerr=pstd.(mid[date][:S]),fmt="+",color="black")
        ax.xaxis.set_major_formatter(majorformatter)
        xlabel("Time")
        ylabel(L"Cross-sectional area $S$ [$\mathrm{m^2}$]")

        ## linear fit and scattering of measurements
        xdata = mid[date][:t_inj]
        xdata = (xdata.-xdata[1])*1e-3./Dates.Millisecond(1)
        ydata = mid[date][:S]

        @. fct(xx,p) = p[1]+p[2]*xx
        para_weights = [0.5,0.5]
        dS_fit = curve_fit_MCMeasurements(fct, xdata, ydata, para_weights) # determine linear trend
        trend1 = pmean(dS_fit[1]).+pmean(dS_fit[2])*xdata
        trend2 = pmean(dS_fit[1])+pstd(dS_fit[1]).+(pmean(dS_fit[2])+pstd.(dS_fit[2]))*xdata
        trend3 = pmean(dS_fit[1])-pstd(dS_fit[1]).+(pmean(dS_fit[2])-pstd.(dS_fit[2]))*xdata
        pl2 = fill_between(mid[date][:t_inj], trend2, trend3, color="lightblue")

        ## S from computed dS_dt
        t_plot = mid[date][:t_inj]
        S_mean = trend1[midpt]
        dt = (t_plot-mid[date][:t_inj][midpt])*1e-3/Dates.Millisecond(1)

        dS_dt = model_runs[date][:dSdt_ct] - model_runs[date][:closure]
        dS_dt_mean = mean(dS_dt[isnan.(pmean.(dS_dt)).==0],dims=1)
        S_calc = S_mean .+ dS_dt_mean .* dt
        plot(mid[date][:t_inj], pmean.(S_calc),"-",color="black")
        pl3 = fill_between(mid[date][:t_inj], pmean.(S_calc)-pstd.(S_calc), pmean.(S_calc)+pstd.(S_calc), color="lightgrey", label="Scatter of model")

        title(date[1:2] * "-Aug")

        if date == "0908"
            legend((custom_lines[1], pl3, pl1, pl2),
                (L"Predicted growth rate, $c_t$-gradient model",
                "Error of prediction",
                "Measurements",
                "Linear fit of measurements"),
                bbox_to_anchor=(0.25, -0.25, 1.5, .102), # position of the legend (x0, y0, width, height)
                ncol = 2,
                )
        end
        text(0.05, 0.95, lab, transform=ax.transAxes, ha="left", va="top")
    end
    gcf()
end

"""
For a chosen parameter of heat transfer, plot values of each experiment on 9 and 21 August to see the variation over the days. All uncertainties with one standard deviation.
Possible parameters: z_eq, tau_eq, tau_w; for tau_w it is possible to add measurements as a third argument
"""
function plot_heat_params_timeresolved(parameters, y_labels, measurements)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 13

    panellabs = [L"\bf{a}", L"\bf{b}", L"\bf{c}", L"\bf{d}", L"\bf{e}", L"\bf{f}", L"\bf{g}", L"\bf{h}", L"\bf{i}", L"\bf{j}", L"\bf{k}", L"\bf{l}", L"\bf{m}"]

    figure(figsize=(10, 12), dpi=pixel_resolution)
    for (np, (parameter, measurement, y_label)) in enumerate(zip(parameters, measurements, y_labels))
        for (nd, (date, tit)) in enumerate(zip(["0908", "2108"], ["9-Aug/AM15", "21-Aug/AM13"]))
            ax = subplot(length(parameters), 2, 2*(np-1)+nd)
            if measurement !== nothing
                fill_between(1:length(parameter[date][:standard]), pquantile(only(measurement[date]), 0.025), pquantile(only(measurement[date]), 0.975), label="Measured", color="grey", alpha=0.3)
            end
            for (corr, col) in zip(keys(parameter[date]), ["green", "pink", "blue", "darkorange", "blueviolet"])
                errorbar(1:length(parameter[date][corr]), pmean.(parameter[date][corr]), yerr=pstd.(parameter[date][corr]), fmt="_", color=col, label=uppercasefirst(string(corr)), markersize=7)
            end
            xticks(2:2:length(parameter[date][:standard]))
            if date == "0908"
                ylabel(y_label)
            end
            if np == 1
                title(tit)
            elseif np == length(parameters)
                xlabel("Index of tracer experiment")
            end
            if date == "0908" && np == length(parameters) # measurement !== nothing
                legend(ncol=3, loc="lower right", framealpha = 0.4,
                   columnspacing=1.0, fontsize=11, handlelength=1.0)
            end
            text(0.0, 1.05, panellabs[2*(np-1)+nd], transform=ax.transAxes, ha="right", va="bottom")
        end
    end
    gcf()
end


"""
Plot the distributions of paper figure 5d)
"""
function plot_tauw_hist(tauw_calc, tauw_meas)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 12

    # custom legend for Patches
    # get the rgb colors with PyPlot.matplotlib.colors.to_rgb("green")
    linew = 1.2
    Patch = PyPlot.matplotlib.patches.Patch
    custom_legend = (Patch(facecolor=(0, 0, 0, 0.2), edgecolor=(0, 0, 0, 1), lw=linew), # black
                     Patch(facecolor=(0, 0.5019607843137255, 0, 0.3), edgecolor=(0, 0.5019607843137255, 0, 1), lw=linew), # green
                     Patch(facecolor=(0.5411764705882353, 0.16862745098039217, 0.8862745098039215, 0.3), edgecolor=(0.5411764705882353, 0.16862745098039217, 0.8862745098039215, 1), lw=linew), # blueviolet
                     Patch(facecolor=(0, 0, 1, 0.3), edgecolor=(0, 0, 1, 1), lw=linew), # blue
                     Patch(facecolor=(1.0, 0.5490196078431373, 0.0, 0.3), edgecolor=(1.0, 0.5490196078431373, 0.0, 1), lw=linew), # darkorange
                     Patch(facecolor=(1.0, 0.7529411764705882, 0.796078431372549, 0.6), edgecolor=(1.0, 0.7529411764705882, 0.796078431372549, 1), lw=linew) # pink
    )

    figure(figsize=(16,7), dpi=pixel_resolution)
    for (nd, (date, lab)) in enumerate(zip(["0908", "2108"], [L"\bf{a}", L"\bf{a}"]))
        ax = subplot(2,1,nd)
        hist(Array(tauw_meas[date]), bins=100, color="black", alpha=0.2, density=true)
        hist(Array(tauw_meas[date]), bins=100, color="black", histtype="step", lw=linew, density=true)
        tick_params("y", left=false, labelleft=false, right=false, labelright=false, direction="in", color="white")
        ylabel("Frequency")
        for (corr, col) in zip(keys(tauw_calc[date]), ["green", "blueviolet", "blue", "darkorange", "pink"])
            ax.twinx()
            hist(Array(mean(tauw_calc[date][corr], dims=1)), bins=100, color=col, alpha=0.3, density=true, zorder=1)
            hist(Array(mean(tauw_calc[date][corr], dims=1)), bins=100, color=col, histtype="step", lw=linew, density=true, zorder=1)
            tick_params("y", left=false, labelleft=false, right=false, labelright=false)
        end
        legend(custom_legend,
               ("Measured", "Standard", "Lunardini", "Vincent", "Ogier", "Gnielinski"),
               loc = "lower left",
               ncol = 2,
               handlelength = 1.5,
               handletextpad = 0.6,
               columnspacing = 1.0
        )
        #text(0.03, 0.9, lab, transform=ax.transAxes, ha="left", va="top")
        text(0.95, 0.9, date[1:2] * "-Aug", transform=ax.transAxes, ha="right", va="top")
        xlim((-0.12, 0.22))
    end
    gcf()
end

"""
Plot the daily mean values of all heat transfer parameters (Nu, z_eq, tau_eq, tau_w, tau_w - tau_eq, T_surf) for each Nu-parameterisation
Fig. 5 in the paper and Fig. 2 in supplements (for the latter use additional argument)
Additional output: text to create latex table of the plotted values
"""
function plot_heat_transfer_params(Nus, z_eq, tau_eq, tau_w, tau_diff, tauw_measured; T_surf=nothing)
    # font properties
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.size"] = 14

    # customise legend
    Patch = PyPlot.matplotlib.patches.Patch
    custom_legend = [Patch(facecolor=(0, 0, 0, 0.2), edgecolor=(0, 0, 0, 1)),
                     Patch(facecolor=(0, 0.5019607843137255, 0, 0.1), edgecolor=(0, 0.5019607843137255, 0, 1))]

    # define Tuples and Arrays to loop over
    paras = (Nu = (L"Nu", Nus),
             z0 = (L"z_{eq}\,(\mathrm{m})", z_eq),
             tau_eq = (L"\tau_{eq}\,(\mathrm{\degree C})", tau_eq),
             tau_w_off = (L"\tau_w\,(\mathrm{\degree C})", tau_w),
             tau_diff = (L"\tau_w-\tau_{eq}\,(\mathrm{\degree C})", tau_diff)
             )
    xticklabs = [:standard, :sommers, :lunardini, :vincent, :ogier, :gnielinski]
    panellabs = [L"\bf{a}", L"\bf{b}", L"\bf{c}", L"\bf{d}", L"\bf{e}", L"\bf{f}"]
    days = ["0908", "2108"]

    # how many standard deviations to plot
    nr_of_stds = 1

    # draw figure
    if T_surf !== nothing
        paras = (; paras... , tau_surf = (L"\tau_{surf}\,(\mathrm{\degree C})", T_surf)) # if T_surf is plotted as well (supplements)
        figure(figsize=(6,12), dpi=pixel_resolution)
    else
        figure(figsize=(6,10), dpi=pixel_resolution)
    end

    table = ""
    for (sp,para) in enumerate(paras)
        out = []
        lab,p = para
        xs = []
        ys = []
        yerr = []

        ax = (nothing for i=1:100)
        for (dx, date) in zip([-0.05,0.05], days) # loop over days
            ax = subplot(length(paras),1,sp)
            o = lab # for latex table
            if date == "2108"
                o = ""
            end
            for (x,coef) in enumerate(xticklabs) # loop over Nu parameterisations
                mean_ = sum(p[date][coef])/length(p[date][coef])
                std_ = max(pstd(mean_), std(pmean.(p[date][coef])))
                mean_ = pmean(mean_)
                if date == "0908"
                    o = o * " & " * string(round(mean_, digits=3) ± round(std_, digits=3))
                elseif date == "2108"
                    o = o * " & \\textbf{" * string(round(mean_, digits=3) ± round(std_, digits=3)) * "}"
                end
                push!(xs, x+dx)
                push!(ys, mean_)
                push!(yerr, nr_of_stds*std_)
            end
            push!(out, o)
        end
        errorbar(xs[1:length(xticklabs)], ys[1:length(xticklabs)], yerr=yerr[1:length(xticklabs)], fmt="k+", markersize=7) # 9 August
        errorbar(xs[length(xticklabs)+1:end], ys[length(xticklabs)+1:end], yerr=yerr[length(xticklabs)+1:end], fmt="g+", markersize=7) # 21 August
        ylabel(lab)
        xticks(1:length(xticklabs), ["" for i=1:length(xticklabs)])
        if p == Nus
            legend(["9-Aug/AM15", "21-Aug/AM13"],
                           loc="upper left",
                           ncol = 1,
                           handletextpad = 0.2,
                           borderpad = 0.2,
                           columnspacing = 0.1,
                           handlelength = 0.7,
                           fontsize = 11
                           )
            yticks([0, 5000, 10000])
        elseif  p == tau_w
            for (col, alph, date) in zip(["black","green"], [0.2, 0.1], days)
                tw = Array(tauw_measured[date][1])
                fill_between([xs[1]-0.06, xs[end]+0.06], [quantile(tw, 0.025), quantile(tw, 0.025)],
                                                         [quantile(tw, 0.975), quantile(tw, 0.975)],
                                                         color=col,
                                                         alpha=alph)
                plot([xs[1]-0.06, xs[end]+0.06], [quantile(tw, 0.975), quantile(tw, 0.975)], color=col, linestyle="-", linewidth=1.0)
                plot([xs[1]-0.06, xs[end]+0.06], [quantile(tw, 0.025), quantile(tw, 0.025)], color=col, linestyle="-", linewidth=1.0)
                legend(custom_legend, ("Measured 9-Aug/AM15", "Measured 21-Aug/AM13"),
                       loc="upper right",
                       ncol = 1,
                       handletextpad = 0.5,
                       borderpad = 0.2,
                       columnspacing = 0.1,
                       handlelength = 1.0,
                       fontsize = 11
                )
            end
        elseif p == tau_diff
            hlines(0, xlim()[1], xlim()[2], "k", linestyle="--", linewidth=0.6)
        elseif p == T_surf
            ylim([-3.5, 1.5])
        end
        text(1.05, 1.0, panellabs[sp], transform=ax.transAxes, ha="left", va="top")
        # latex table
        out = join(out, "\\\\ \n")
        table = table * out * "\\\\ \n"
    end
    xticks(1:length(xticklabs), [uppercasefirst(string(xl)) for xl in xticklabs], rotation = -40, rotation_mode="anchor", ha="left")
    tight_layout()
    return gcf(), table
end

# Functions or making map (Fig. 1 in paper)
"Get map extent (https://github.com/JuliaPlots/Plots.jl/issues/394)"
get_map_extent(p) = round.(Int, (p[1].o.get_xlim())), round.(Int, (p[1].o.get_ylim()))
"Get map center (https://github.com/JuliaPlots/Plots.jl/issues/394)"
get_map_center(p) = round.(Int, (mean(p[1].o.get_xlim()), mean(p[1].o.get_ylim())) )

end # module
