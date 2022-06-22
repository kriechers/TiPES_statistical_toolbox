
################################################################################
# this program identifies the midpoint of rapid transitions during the last
# glacial period (D/O events) of your record of interest
#
# the identified transitions are related to NGRIP detected D/O events
# (Rasmussen et al., 2014)
#
# The identified D/O events are then stacked and the mean offset of stadial and
# interstadial conditions is provided as output.
#
# This code is part of the TiPES deliverable 1.2 (Statistacal toolbox for
# the detection of tipping points
################################################################################

# using Pkg
#
# Pkg.add("Plots")
# Pkg.add("DelimitedFiles")
# Pkg.add("Statistics")
# Pkg.add("Interpolations")
# Pkg.add("Infiltrator")      # shall enable a matlab like 'keyboard' comand: @infiltrate

using Plots, Interpolations
using DelimitedFiles, Statistics


################################################################################
### two neighboring 100 a windows, their runnung mean and the resulting difference
################################################################################
function running_mean_diff(time_series, age, window)

    itp = LinearInterpolation(age,time_series)
    age_interp = age[1]:1:age[end]
    interp_dummy = itp(age_interp)
    age_interp =age_interp[findall(x->x==false, isnan.(interp_dummy))]
    interp_dummy =interp_dummy[findall(x->x==false, isnan.(interp_dummy))]
    differenz = zeros(length(interp_dummy))

    if length(interp_dummy) < 2*window
        window = floor(Int,length(interp_dummy)/2)
    end
    for i = 11:window
        mean_pre = mean(interp_dummy[1:i-1])
        mean_post = mean(interp_dummy[i:i+window])
        differenz[i] = mean_pre-mean_post # positive values show increase of d18O
    end

    for i = 1+window : length(interp_dummy)-window
        mean_pre = mean(interp_dummy[i-window:i-1])
        mean_post = mean(interp_dummy[i:i+window])
        differenz[i] = mean_pre-mean_post # positive values show increase of d18O
    end

    for i = length(interp_dummy)-window+1 : length(interp_dummy)-10
        mean_pre = mean(interp_dummy[i-window:i])
        mean_post = mean(interp_dummy[i:end])
        differenz[i] = mean_pre-mean_post # positive values show increase of d18O
    end

    return differenz, interp_dummy, age_interp

end


############################################################################
### reading data files
cd(@__DIR__) # to change to path of this *.jl file
proxyfile = "YourRecord.txt"    # for rows:   depth  age  age error  proxy (tab separated)
agefile = "YourDatingTable.txt" # three rows: depth  age  age error (tab separated)
record_age = readdlm(proxyfile, Float64, comments = true, comment_char = '#')
DO_direction = 1 # 1/-1 : proxy values for stad lower than interstad/opposite
age = readdlm(agefile, Float64, comments = true, comment_char = '#')

NGRIP = readdlm("NGRIP.txt", header = true, '\t', Float64, '\n')
NGRIP_events = readdlm("NGRIP_StartOfGImainmain.txt", '\t', '\n')
############################################################################


### extract proxy and, age data ########################################
age_dummy = record_age[:,2]    # removing those wrongly inserted data in SISAL
ageerr_dummy = record_age[:,3]  # removing those wrongly inserted data in SISAL
d18O_dummy = record_age[:,4]    # removing those wrongly inserted data in SISAL
############################################################################

### sort in time-increasing order ##########################################
p = sortperm(age_dummy)
d18O_dummy = d18O_dummy[p]
ageerr_dummy = ageerr_dummy[p]
age_dummy = age_dummy[p]
############################################################################

### for figures: fix x range for stal and NGRIP (on the nearest millennium)
x_low = floor(minimum(age_dummy)./1000).*1000
x_high = ceil(maximum(age_dummy)./1000).*1000
############################################################################

### take NGRIP Stadial-Interstadial transitions (from Rasmussen et al., 2014)
### one by one and look for appropriate transitions in individual spelothems
if maximum(age_dummy) > maximum(NGRIP_events[:,2])
    id_DO_high = length(NGRIP_events[:,2])
else
    id_DO_high = findlast(x -> x < maximum(age_dummy), NGRIP_events[:,2])
end
if minimum(age_dummy) < minimum(NGRIP_events[:,2])
    id_DO_low = 1
else
    id_DO_low = findfirst(x -> x > minimum(age_dummy), NGRIP_events[:,2])
end

global timing_DO = zeros(0)
global p3 = plot(size=[600,240], seriestype=:line,
    yguidefont = font(:blue),
    ytickfont = font(:blue))

for i = id_DO_low:id_DO_high
    if i == 1
        tau=abs(NGRIP_events[i,2]-NGRIP_events[i+1,2])
    elseif i == length(NGRIP_events[:,2])
        tau = abs(NGRIP_events[i,2]-NGRIP_events[i-1,2])
    else
        tau = minimum([abs(NGRIP_events[i,2]-NGRIP_events[i-1,2]),abs(NGRIP_events[i,2]-NGRIP_events[i+1,2])])
    end

    ###find age errors from U/Th ages in dating.csv
    if age == empty
        println("do something here")
    else
        dummy_age = age[2]#dating[findall(x -> x == sampleDO.entity_id[i2[m]], dating.entity_id),:corr_age]
        dummy_ageerr = age[3]#dating[findall(x -> x == sampleDO.entity_id[i2[m]], dating.entity_id),:corr_age_uncert_pos]
        id_dummy_age = findall(x -> x == minimum(abs.(NGRIP_events[i,2].-dummy_age)),abs.(NGRIP_events[i,2].-dummy_age))
    end
    ### search in between +/-tau/2 +/-age_err of stal around each D/O event
    low = NGRIP_events[i,2] .- tau/2 .- dummy_ageerr[id_dummy_age[1]]
    high = NGRIP_events[i,2] .+ tau/2 .+ dummy_ageerr[id_dummy_age[1]]
    d18O_dummy1= d18O_dummy[low .< age_dummy .< high]
    age_dummy1 = age_dummy[low .< age_dummy .< high]

    ############################################################################
    ### two 'ii' year long running means and their difference for d18O --> differenz
    ### do some detrending before, if necessary
    ############################################################################
    if isempty(age_dummy1) || length(age_dummy1) == 1 # this is necessary in case of a growth stop
        #println(i)
    else
        global ii = minimum([200,maximum([floor(Int,(age_dummy1[end]-age_dummy1[1])/2/100)*100,50])])
        while ii >= 50
            global (differenz, d18O_interp_dummy, age_interp) = running_mean_diff(d18O_dummy1.*DO_direction[1],age_dummy1,200)
            d18O_interp_dummy = d18O_interp_dummy.*DO_direction[1]
                #'.*DO_direction[1][m,3]' to account for direction of DO response in d18O

            ### find local minima in differenz above (below) std ######
            global streuung_diff_AR1_O = 1*std(differenz[findall(x->x.!=0,diff(differenz))]) # to neglect influence of hiatus
            if isnan(streuung_diff_AR1_O) global streuung_diff_AR1_O = 0 end
            jj = findall(x -> x > streuung_diff_AR1_O, differenz)
            global jj1 = findall(x -> x != 1, diff(jj))

            if isempty(jj1) && isempty(jj)
                global timing_dummy = zeros(1)
                timing_dummy_length = zeros(1)
                max_DO_offset = zeros(1)
            elseif isempty(jj1)
                global timing_dummy = zeros(1)
                timing_dummy_length = zeros(1)
                max_DO_offset = zeros(1)
                timing_dummy[1] = age_interp[jj[1] .+
                    findall(x -> x .== maximum(differenz[jj[1]:jj[end]]),differenz[jj[1]:jj[end]])][1]
                timing_dummy_length[1] = age_interp[jj[end]]-age_interp[jj[1]]
                max_DO_offset[1] = maximum(differenz[jj[1]:jj[end]]) - streuung_diff_AR1_O

            else
                global timing_dummy = zeros(length(jj1)+1)
                timing_dummy_length = zeros(length(jj1)+1)
                max_DO_offset = zeros(length(jj1)+1)
                timing_dummy[1] = age_interp[findall(x -> x .== maximum(differenz[1:jj[jj1[1]]]),differenz[1:jj[jj1[1]]])][1]
                timing_dummy_length[1] = age_interp[jj[jj1[1]]]-age_interp[jj[1]]
                max_DO_offset[1] = maximum(differenz[1:jj[jj1[1]]]) -
                    streuung_diff_AR1_O
                for i = 2:length(jj1)
                    timing_dummy[i] = age_interp[jj[jj1[i-1]+1]-1 +
                        findall(x -> x .== maximum(differenz[jj[jj1[i-1]+1]:jj[jj1[i]]]),differenz[jj[jj1[i-1]+1]:jj[jj1[i]]])[1]]
                    timing_dummy_length[i] = age_interp[jj[jj1[i]]]-age_interp[jj[jj1[i-1]+1]]
                    max_DO_offset[i] = maximum(differenz[jj[jj1[i-1]+1]:jj[jj1[i]]]) - streuung_diff_AR1_O
                end
                timing_dummy[end] = age_interp[jj[jj1[end]+1]-1 +
                    findall(x -> x .== maximum(differenz[jj[jj1[end]+1]:end]),differenz[jj[jj1[end]+1]:end])[1]]
                timing_dummy_length[end] = age_interp[jj[end]]-age_interp[jj[jj1[end]+1]]
                max_DO_offset[end] = maximum(differenz[jj[jj1[end]+1]:end]) - streuung_diff_AR1_O

            end

            if mean(timing_dummy_length)<=500
                timing_dummy = timing_dummy[max_DO_offset.==maximum(max_DO_offset)]
                if isempty(timing_dummy)
                    timing_dummy = zeros(1)
                end
                append!(timing_DO, timing_dummy)
                break
            end
            if ii == 50
                break
            end

            global ii = ii-50
        end

        ### plot differenz and plot found local minima #####################
        global p3 = plot!(age_interp,differenz,linecolor=:blue,
        xlabel = "time [ka BP]", ylabel= "proxy [unit]",legend = false,
        title = "difference of two neighboring running windows",
        titlefontsize = 12, xlims = (x_low,x_high),xflip = true,framestyle = :box,
        xguidefont = font(10),xtickfont = font(8),
        yguidefont = font(:blue,10),ytickfont = font(:blue,8))

        p3 = plot!([age_dummy1[1],age_dummy1[end]],[streuung_diff_AR1_O,streuung_diff_AR1_O],linecolor=:cyan)
        p3 = vline!([convert(Array{Float64},(NGRIP_events[:,2]))], linecolor = :green)
        ####################################################################
    end
end
############################################################################

### remove double non-zero entries ('double': defined as smaller distance as 200a)
# distance matrix dist
println(timing_DO)
for dummyyy = 1:1
    z = 0
    while z == 0
        dist = zeros(length(timing_DO)-1)
        for i = 1:length(timing_DO)-1
            dist[i] = abs(timing_DO[i+1]-timing_DO[i])
        end

        hh = findall(x-> 0 .< x .< 200, dist)
        if isempty(hh)
            z=1
        else
            for j1 = 1:length(hh)
                timing_DO[hh[j1][1]] = timing_DO[hh[j1][2]]
            end
        end

        unique!(timing_DO)
        filter!(x->x!=0,timing_DO)
    end
end
if isempty(timing_DO) else p3 = vline!(timing_DO, linecolor = :blue) end
##############################################################


### plot proxy (blue) ########################################
p2 = plot( age_dummy, d18O_dummy,
    legend = false,
    titlefontsize = 12, size=[600,240], seriestype=:line,
    xlabel = "time [ka BP]", ylabel= "proxy [unit]",
    xflip = true, linecolor = :blue, yguidefont = font(:blue,10),
    ytickfont = font(:blue,8),framestyle = :box,
    xguidefont = font(10),xtickfont = font(8), xlims = (x_low,x_high))
if isempty(timing_DO) else p2 = vline!(timing_DO, linecolor = :blue) end
############################################################################

### plot NGRIP data and the defined starting point of D/O events (green) ###
p1 = plot(NGRIP[1][:,1].-50, NGRIP[1][:,2],framestyle = :box,
xlabel = "time [ka BP]", ylabel= "d18O [o/oo VSMOW]",legend = false,
title = "NGRIP ice core data", titlefontsize = 12, ymirror = true,
xlims = (x_low,x_high), xflip = true, color = :lightgreen,
yguidefont = font(:lightgreen,10),ytickfont = font(:lightgreen,8),
xguidefont = font(10),xtickfont = font(8))
p1 = vline!([convert(Array{Float64},(NGRIP_events[:,2]))], linecolor = :green)
############################################################################

### take the transition and +/- 300 years, then put all on the same d18O-level
### before the event and stack
if isempty(timing_DO)
    p4=plot(300:-1:-300, [300:-1:-300,-300:1:300],
    xlabel = "relative time [years]", ylabel= "d18O [permil VPDB]",
    framestyle = :box, tickfont = font(6), guidefontsize = 6,
    titlefontsize = 6, legend = false, color = :blue)
    annotate!(-200,-250,text("No DOs detected", color = :blue, :left, 8))
else
    x_dummy = zeros(601,length(timing_DO))
    save_d18O = zeros(601,length(timing_DO))
    for i = 1:length(timing_DO)
        period = round.(Int,age_dummy[findall(x -> -500<=x.-timing_DO[i]<=500,age_dummy)]).-round(Int,timing_DO[i])
        x_dummy[1:length(period),i] = -period #'-'sign to ensure that everything smaller than 0 is stadial
        save_d18O[1:length(period),i] = d18O_dummy[findall(x -> -500<=x.-timing_DO[i]<=500,age_dummy)]
    end

    global last=findlast(x -> x !=0, x_dummy[:,1])
    global p4 = plot(x_dummy[1:last,1], save_d18O[1:last,1], color = :grey,
        xlabel = "time relative to transition [years]", ylabel= "proxy [unit]",
        framestyle = :box, tickfont = font(6), guidefontsize = 6,
        titlefontsize = 6,xlims=(-300,300))

    for i = 2:length(timing_DO)
        last=findlast(x -> x !=0, x_dummy[:,i])
        p4 = plot!(x_dummy[1:last,i], save_d18O[1:last,i], legend = false, color = :grey)
    end
end

plot(p4,size=[300,200],clearfig=false)
savefig("TimeSeries_Transitions.pdf")


# offset: interstadial - stadial
off18O = zeros(length(timing_DO))   #mean over 250a
for i = 1:length(timing_DO)
    off18O[i] = mean(save_d18O[(x_dummy[:,i].!=0) .& (50 .< x_dummy[:,i] .< 300),i]) -
                mean(save_d18O[(x_dummy[:,i].!=0) .& (-300 .< x_dummy[:,i] .< -50),i])
end


# plot offset
p6 = plot(timing_DO, off18O, color = :black, seriestype = :scatter,
    title = "interstadial-stadial difference",
    xlabel = "time [ka]", ylabel= "proxy [Delta unit]", xflip = true,
    framestyle = :box, tickfont = font(8), guidefontsize = 10,
    titlefontsize = 12, legend = false, ymirror = true, xlims = (x_low,x_high))
### plot all on one page and save as pdf ###################################
plot(p2,p1,p3,p6, layout=grid(4,1,heights=[0.4,0.2,0.2,0.2]),size=[600,600],clearfig=false)
savefig("AllTransitions.pdf")
############################################################################


### save data
open("DetectedEvents.txt"; write=true) do f
    write(f, "# time [ka]   IS-S proxy difference [unit]\n")
    writedlm(f, [round.(timing_DO) round.(off18O,digits=2)])
end
