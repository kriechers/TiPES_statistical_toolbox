
################################################################################
# this program identifies D/O events of some predefined records from SISAL2
# The identified D/O events are then stacked and the mean offset of stadial and
# interstadial conditions is provided as output.
#
# This work is intended to be published in ..., providing a blueprint for
# climate models
################################################################################

# Comment blocks: 'ctrl' + '/'

# using Pkg
#
# Pkg.add("CSV")
# Pkg.add("DataFrames")
# Pkg.add("DataFramesMeta")
# Pkg.add("Plots")
# Pkg.add("GMT")
# Pkg.add("Interpolations")
# Pkg.add("Downloads")
# Pkg.add("InfoZIP")
# Pkg.add("ZipFile")
# Pkg.add("Infiltrator")      # shall enable a matlab like 'keyboard' comand: @infiltrate

using CSV, Plots, Interpolations
using DataFrames, DelimitedFiles, Statistics
using ZipFile, InfoZIP, Downloads


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

################################################################################
### find continent of cave location
################################################################################
function find_continent(lat, lon)

    x = zeros(length(lat))
    for i = 1:length(lat)
        # Europe
        if 35 <= lat[i] <= 65 && -10 <= lon[i] <= 40
            x[i] = 1
        # North America
        elseif 10 <= lat[i] <= 60 && -140 <= lon[i] <= -55
            x[i] = 2
        # South America
        elseif -50 <= lat[i] <= 10 && -80 <= lon[i] <= -35
            x[i] = 3
        # Africa
        elseif -35 <= lat[i] <= 35 && -15 <= lon[i] <= 60
            x[i] = 4
        # China, India (Asia?)
        elseif 15 <= lat[i] <= 60 && 65 <= lon[i] <= 125
            x[i] = 5
        # Oceania (including Australia)
        elseif -50 <= lat[i] <= 15 && 95 <= lon[i] <= 180
            x[i] = 6
        # other parts of the world
        else
            x[i] = 7
        end
    end

    return x

end # function
################################################################################

cd(@__DIR__) # to change to path of this *.jl file

### download SISALv2.zip and unpack into "./SISALv2_csv/"
if !isdir("./SISALv2_csv")
    mkdir("./SISALv2_csv/")
    Downloads.download("https://researchdata.reading.ac.uk/256/31/sisalv2_csv.zip", "sisalv2_csv.zip")
    InfoZIP.unzip("sisalv2_csv.zip", "./SISALv2_csv/")
    rm("sisalv2_csv.zip")
end

### reading data files
cd("SISALv2_csv") # to change the path relative to where it was before
entity = CSV.read("entity.csv", DataFrame, missingstring = "NA")
gap = CSV.read("gap.csv", DataFrame, types=(Int64,String)) # reads CSV file relative to where the julia is running
composite_link_entity = CSV.read("composite_link_entity.csv", DataFrame)
d13C = CSV.read("d13C.csv", DataFrame, missingstring = "NA")
d18O = CSV.read("d18O.csv", DataFrame, missingstring = "NA")
dating_lamina = CSV.read("dating_lamina.csv", DataFrame)
dating = CSV.read("dating.csv", DataFrame, missingstring = "NA")#treating all NULL as missing values
rename!(dating,Symbol.([:dating_id,
    :entity_id,:date_type,:depth_dating,:dating_thickness,:lab_num,:material_dated,
    :min_weight,:max_weight,:uncorr_age,:uncorr_age_uncert_pos,:uncorr_age_uncert_neg,
    :C14_correction,:calib_used,:date_used,:c238U_content,:c238U_uncertainty,
    :c232Th_content,:c232Th_uncertainty,:c230Th_content,:c230Th_uncertainty,
    :a230Th_232Th_ratio,:a230Th_232Th_ratio_uncertainty,:a230Th_238U_activity,
    :a230Th_238U_activity_uncertainty,:a234U_238U_activity,:a234U_238U_activity_uncertainty,
    :ini_230Th_232Th_ratio,:ini_230Th_232Th_ratio_uncertainty,:decay_constant,
    :corr_age,:corr_age_uncert_pos,:corr_age_uncert_neg,:date_used_lin_interp,
    :date_used_lin_reg, :date_used_Bchron,:date_used_Bacon,:date_used_OxCal,
    :date_used_copRa, :date_used_StalAge]))    # it is necessary to rename those lines with a number on first position
entity_link_reference = CSV.read("entity_link_reference.csv", DataFrame)
hiatus = CSV.read("hiatus.csv", DataFrame)
notes = CSV.read("notes.csv", DataFrame)
original_chronology = CSV.read("original_chronology.csv", DataFrame)
reference = CSV.read("reference.csv", DataFrame)
sample = CSV.read("sample.csv", DataFrame, missingstring = "NA",copycols=true)
sisal_chronology = CSV.read("sisal_chronology.csv", DataFrame)
site = CSV.read("site.csv", DataFrame)
cd(@__DIR__)

NGRIP = readdlm("NGRIP.txt", header = true, '\t', Float64, '\n')
DO_direction = readdlm("DO_direction.txt", header = true, '\t', '\n')
NGRIP_events = readdlm("NGRIP_StartOfGImainmain.txt", '\t', '\n')

################################################################################
#### find all speleothems as defined in DO_direction
################################################################################
#sampleDO = sample[findall(in(DO_direction[1][:,1]),sample.entity_id),:]
sampleDO = DataFrame()
for i = 1:length(DO_direction[1][:,1])
    append!(sampleDO, sample[findall(x->x==DO_direction[1][i,1],sample.entity_id),:])
end

### all stalagmites in sampleDO
stalDO = sort(unique(sampleDO,:entity_id),:entity_id)
println(size(stalDO,1)," stalagmites avaialable in SISAL covering the requested period.")

### all cave sites with stalagmites with DOs
#siteDO = entity[findall(in(stalDO[!,:entity_id]),entity.entity_id),:]
siteDO = DataFrame()
for i = 1:size(stalDO)[1]
    append!(siteDO, entity[findall(x->x.==stalDO[i,:entity_id],entity.entity_id),:])
end
siteDO = sort(unique(siteDO,:site_id),:site_id)
println("Those stalagmites are from $(size(siteDO,1)) caves around the world.")

# find different stalagmites of sampleDO
i = findall(x->x!=0,diff(sampleDO.entity_id)).+1
i1= zeros(length(i)+2)
i1[1]=1
i1[2:length(i)+1] = i
i1[length(i)+2] = length(sampleDO.entity_id)+1
i1=convert(Vector{Int},i1)      # i1 provides first index of a stal (except of i1(length(i1)) )

sort!(sampleDO,[:entity_id, :sample_id])
println(size(sampleDO))
sampleDO = dropmissing(sampleDO, :mineralogy)
println(size(sampleDO))

# find different stalagmites within sampleDO
i = findall(x->x!=0,diff(sampleDO.entity_id)).+1
i2= zeros(length(i)+2)
i2[1]=1
i2[2:length(i)+1] = i
i2[length(i)+2] = length(sampleDO.entity_id)+1
i2=convert(Vector{Int},i2)      # i2 provides first index of a stal

mean_off18O = zeros(length(i2)-1)
std_off18O = zeros(length(i2)-1)
df = DataFrame(entity_id = Int[], entity_name = String[], site_id = Int[], site_name = String[],
    latitude = Float64[], longitude = Float64[], distance_NGRIP = Float64[],
    timing = Int[], mean_off18O = Float64[], continent = Float64[])
radius_earth = 6371
coord_NGRIP = [75.1, -42.32, radius_earth] # lat,lon in degree
coord_NGRIP = [radius_earth*cosd(coord_NGRIP[1])*cosd(coord_NGRIP[2]),
               radius_earth*cosd(coord_NGRIP[1])*sind(coord_NGRIP[2]),
               radius_earth*sind(coord_NGRIP[1])] # in coordinates

distance_NGRIP = zeros(length(i2)-1)
xxx = zeros(length(i2)-1)

for m = 1:length(i2)-1

    ### info for plot title ####################################################
    idx = findall(x -> x == sampleDO.entity_id[i2[m]],entity.entity_id)
    idx1 = findall(x -> x == entity.site_id[idx[1]],site.site_id)
    sample_id_dummy = sampleDO.sample_id[i2[m]:i2[m+1]-1]
    ############################################################################
    println(m," ",entity.entity_name[idx][1])


    ### extract d18O and age data ########################################
    d18O_dummy = zeros(length(sample_id_dummy))
    age_dummy = zeros(length(sample_id_dummy))

    for n = 1:length(sample_id_dummy)
        dummy1 = findall(x-> x.==sample_id_dummy[n],d18O.sample_id)
        if isempty(dummy1)
        else
            d18O_dummy[n] = d18O[d18O.sample_id .== sample_id_dummy[n],:d18O_measurement][1]
        end
        dummy = findall(x-> x.==sample_id_dummy[n],original_chronology.sample_id)
        if isempty(dummy)
        else
            age_dummy[n] = original_chronology[original_chronology.sample_id .== sample_id_dummy[n],:interp_age][1]
        end
    end
    age_dummy = age_dummy[d18O_dummy.!=0.0]    # removing those wrongly inserted data in SISAL
    d18O_dummy = d18O_dummy[d18O_dummy.!=0.0]    # removing those wrongly inserted data in SISAL

    ############################################################################

    ### sort in time-increasing order ##########################################
    p = sortperm(age_dummy)
    d18O_dummy = d18O_dummy[p]
    age_dummy = age_dummy[p]
    d18O_dummy=d18O_dummy[20000 .< age_dummy .< 120000]
    age_dummy=age_dummy[20000 .< age_dummy .< 120000].+0.0001
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
        dummy_age = dating[findall(x -> x == sampleDO.entity_id[i2[m]], dating.entity_id),:corr_age]
        dummy_ageerr = dating[findall(x -> x == sampleDO.entity_id[i2[m]], dating.entity_id),:corr_age_uncert_pos]
        dummy_age = dummy_age[findall(!ismissing,dummy_age)] #remove all missing values
        dummy_ageerr = dummy_ageerr[findall(!ismissing,dummy_ageerr)] #remove all missing values
        if isempty(dummy_age) dummy_age=zeros(1) end # in case of no U/Th dates due to composite
        if isempty(dummy_ageerr) dummy_ageerr = 500 end # in case of no U/Th dates due to composite
        id_dummy_age = findall(x -> x == minimum(abs.(NGRIP_events[i,2].-dummy_age)),abs.(NGRIP_events[i,2].-dummy_age))

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
                global (differenz, d18O_interp_dummy, age_interp) = running_mean_diff(d18O_dummy1.*DO_direction[1][m,3],age_dummy1,200)
                d18O_interp_dummy = d18O_interp_dummy.*DO_direction[1][m,3]
                    #'.*DO_direction[1][m,3]' to account for direction of DO response in d18O

                ### find local minima in differenz above (below) 3 times of streuung ######
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
            xlabel = "time [ka BP]", ylabel= "d18O [o/oo VPDB]",legend = false,
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
    z = 0
    while z == 0
        dist = zeros(length(timing_DO),length(timing_DO))
        for i = 1:length(timing_DO)-1
            for j = i : length(timing_DO)
                dist[i,j] = abs(timing_DO[i]-timing_DO[j])
            end
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
        println(timing_DO)
    end
    if isempty(timing_DO) else p3 = vline!(timing_DO, linecolor = :blue) end

    ### plot d18O (blue)  ########################################
    p2 = plot( age_dummy, d18O_dummy,
        legend = false, title = "$(sampleDO.entity_id[i2[m]]), $(entity.entity_name[idx][1]), $(site.site_name[idx1][1])",
        titlefontsize = 12, size=[600,240], seriestype=:line,
        xlabel = "time [ka BP]", ylabel= "d18O [o/oo VPDB]",
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
        title = "SISAL entity_id $(sampleDO.entity_id[i2[m]]), Stalagmite $(entity.entity_name[idx][1]),
        $(site.site_name[idx1][1]) (lat = $(site.latitude[idx1][1]), lon = $(site.longitude[idx1][1]))",
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

        last=findlast(x -> x !=0, x_dummy[:,1])
        p4 = plot(x_dummy[1:last,1], save_d18O[1:last,1], color = :grey,
            title = "SISAL entity_id $(sampleDO.entity_id[i2[m]]), Stalagmite $(entity.entity_name[idx][1]),
            $(site.site_name[idx1][1]) (lat = $(site.latitude[idx1][1]), lon = $(site.longitude[idx1][1]))",
            xlabel = "time relative to transition [years]", ylabel= "d18O [permil VPDB]",
            framestyle = :box, tickfont = font(6), guidefontsize = 6,
            titlefontsize = 6,xlims=(-300,300))

        for i = 2:length(timing_DO)
            last=findlast(x -> x !=0, x_dummy[:,i])
            p4 = plot!(x_dummy[1:last,i], save_d18O[1:last,i], legend = false, color = :grey)
        end

    end

    plot(p4,size=[300,200],clearfig=false)
    savefig("zzz18$m$(entity.entity_name[idx][1]).pdf")


    ############################################################################
    ###statistics

    # offset: interstadial - stadial
    off18O = zeros(length(timing_DO))   #mean over 250a

    for i = 1:length(timing_DO)

        ### offset determination
        off18O[i] = mean(save_d18O[(x_dummy[:,i].!=0) .& (-300 .< x_dummy[:,i] .< -50),i])-
                    mean(save_d18O[(x_dummy[:,i].!=0) .& (50 .< x_dummy[:,i] .< 300),i])
    end

    mean_off18O[m] = mean(off18O)
    std_off18O[m] = std(off18O)

    ############################################################################
    ### calculate distance between NGRIP and Cave location
    coord_cave = [site.latitude[idx1][1],site.longitude[idx1][1],radius_earth]
    coord_cave = [radius_earth*cosd(coord_cave[1])*cosd(coord_cave[2]),
                  radius_earth*cosd(coord_cave[1])*sind(coord_cave[2]),
                  radius_earth*sind(coord_cave[1])]

    angle1 = acos.( (coord_NGRIP[1]*coord_cave[1] .+ coord_NGRIP[2]*coord_cave[2] .+ coord_NGRIP[3]*coord_cave[3]) ./
                 sqrt.(coord_NGRIP[1].^2 .+ coord_NGRIP[2].^2 .+ coord_NGRIP[3].^2) ./
                 sqrt.(coord_cave[1].^2 .+ coord_cave[2].^2 .+ coord_cave[3].^2))
    distance_NGRIP[m] = radius_earth*angle1
    ############################################################################


    ############################################################################
    ### find continents
    ### 1 = Europe, 2 = NA, 3 = SA, 4 = Africa, 5 = Asia, 6 = Oceania/Australia 7 = other
    xx = find_continent(site.latitude[idx1][1], site.longitude[idx1][1])
    xxx[m] = xx[1]
    ############################################################################


    p6 = plot(timing_DO, off18O, color = :black, seriestype = :scatter,
        title = "stadial-interstadial difference",
        xlabel = "time [ka]", ylabel= "d18O [o/oo VPDB]", xflip = true,
        framestyle = :box, tickfont = font(8), guidefontsize = 10,
        titlefontsize = 12, legend = false, ymirror = true, xlims = (x_low,x_high))
    ### plot all on one page and save as pdf ###################################
    plot(p2,p1,p3,p6, layout=grid(4,1,heights=[0.4,0.2,0.2,0.2]),size=[600,600],clearfig=false)
    #plot(p3,p2,p1, layout=grid(3,1,heights=[0.35,0.5,0.15]),size=[600,600],clearfig=false)
    savefig("aaa$(m)$(entity.entity_name[idx][1]).pdf")
    ############################################################################


    ############################################################################

    ############################################################################
    ##save that together with lat,lon,name,SISAL_id etc

    if isempty(timing_DO)

    else
        for i = 1:length(timing_DO)
            push!(df,(sampleDO.entity_id[i2[m]],entity.entity_name[idx][1], site.site_id[idx1][1],
                site.site_name[idx1][1], site.latitude[idx1][1], site.longitude[idx1][1],
                distance_NGRIP[m], round(Int,timing_DO[i]), off18O[i], xxx[m]))
        end
    end
end

################################################################################
### distance-offset -> events
plot(df.distance_NGRIP[df.continent.==1],abs.(df.mean_off18O[df.continent.==1]),color = :black,
    seriestype = :scatter, label = "Europe", markersize = 3)
plot!(df.distance_NGRIP[df.continent.==2],abs.(df.mean_off18O[df.continent.==2]),color = :yellow,
    seriestype = :scatter, label = "N America", markersize = 3)
plot!(df.distance_NGRIP[df.continent.==3],abs.(df.mean_off18O[df.continent.==3]), color = :red,
    seriestype = :scatter, label = "S America", markersize = 3)
plot!(df.distance_NGRIP[df.continent.==4],abs.(df.mean_off18O[df.continent.==4]), color = :orange,
    seriestype = :scatter, label = "Africa", markersize = 3)
plot!(df.distance_NGRIP[df.continent.==5],abs.(df.mean_off18O[df.continent.==5]), color = :blue,
    seriestype = :scatter,label = "Asia", markersize = 3)
display(plot!(df.distance_NGRIP[df.continent.==6],abs.(df.mean_off18O[df.continent.==6]), color = :green,
    seriestype = :scatter, label = "Oceania/Australia", markersize = 3,
    xlabel = "distance from NGRIP [km]", ylabel= "stadial-interstadial d18O difference [o/oo VPDB]",
    framestyle = :box, legend = :topright))
savefig("offset_distance_events.pdf")
savefig("offset_distance_events.png")


################################################################################
### lat-offset -> events
plot(df.mean_off18O[df.continent.==1],df.latitude[df.continent.==1],color = :black,
    seriestype = :scatter, label = "Europe", markersize = 3)
plot!(df.mean_off18O[df.continent.==2],df.latitude[df.continent.==2],color = :yellow,
    seriestype = :scatter, label = "N America", markersize = 3)
plot!(df.mean_off18O[df.continent.==3],df.latitude[df.continent.==3], color = :red,
    seriestype = :scatter, label = "S America", markersize = 3)
plot!(df.mean_off18O[df.continent.==4],df.latitude[df.continent.==4], color = :orange,
    seriestype = :scatter, label = "Africa", markersize = 3)
plot!(df.mean_off18O[df.continent.==5],df.latitude[df.continent.==5], color = :blue,
    seriestype = :scatter,label = "Asia", markersize = 3)
annotate!(3.1,45, text("a)", color = :black, :left, 12))
display(plot!(df.mean_off18O[df.continent.==6],df.latitude[df.continent.==6], color = :green,
    seriestype = :scatter, label = "Oceania/Australia", markersize = 3,
    ylabel = "latitude [°]", xlabel= "stadial-interstadial d18O difference  [o/oo VPDB]",
    framestyle = :box, legend = :bottomright))

savefig("offset_lat_events.pdf")
savefig("offset_lat_events.png")





CSV.write("stack_offset_events.csv",df)
sort!(df,[:site_id,:entity_id])    # order stal after site_id and then after entity_id
        # makes it easier to calculate the number of D/O events- weighted mean

### df2 provides the weighted d18O offset from all speleothems of an inidividual cave
df2 = DataFrame(site_id = Int[], site_name = String[],
    latitude = Float64[], longitude = Float64[], number = Int[], mean_off18O = Float64[],
    std_off18O = Float64[], dist_greenland = Float64[],
    continent = Float64[])
for i = 1:maximum(df.site_id)
    x = df.mean_off18O[df.site_id .== i]
    continent = unique(df.continent[df.site_id .== i])
    distance_greenland = unique(df.distance_NGRIP[df.site_id .== i])
    w = sum(df.site_id .== i)
    s = w
    if sum(w) !== 0
        x_filter = filter(!isnan,x)
        println(i, x, ' ', x_filter, ' ', w)
        #w_filter = filter(x -> x !== 0,w)
        if isempty(x_filter)
            cavemean_d18O = NaN
        else
            cavemean_d18O = median(x_filter)
        end
        cavestd_d18O = std(x_filter)
        if isnan(cavemean_d18O)

        else
            push!(df2,(df.site_id[df.site_id .== i][1], df.site_name[df.site_id .== i][1],
                df.latitude[df.site_id .== i][1],df.longitude[df.site_id .== i][1],
                s,cavemean_d18O,cavestd_d18O,distance_greenland[1], continent[1]))
        end
    end
end

################################################################################
### distance-offset -> caves
plot(df2.dist_greenland[df2.continent.==1],abs.(df2.mean_off18O[df2.continent.==1]),color = :black,
    seriestype = :scatter, label = "Europe", markersize = 3)
plot!(df2.dist_greenland[df2.continent.==2],abs.(df2.mean_off18O[df2.continent.==2]),color = :yellow,
    seriestype = :scatter, label = "N America", markersize = 3)
plot!(df2.dist_greenland[df2.continent.==3],abs.(df2.mean_off18O[df2.continent.==3]), color = :red,
    seriestype = :scatter, label = "S America", markersize = 3)
plot!(df2.dist_greenland[df2.continent.==4],abs.(df2.mean_off18O[df2.continent.==4]), color = :orange,
    seriestype = :scatter, label = "Africa", markersize = 3)
plot!(df2.dist_greenland[df2.continent.==5],abs.(df2.mean_off18O[df2.continent.==5]), color = :blue,
    seriestype = :scatter,label = "Asia", markersize = 3)
display(plot!(df2.dist_greenland[df2.continent.==6],abs.(df2.mean_off18O[df2.continent.==6]), color = :green,
    seriestype = :scatter, label = "Oceania/Australia", markersize = 3,
    xlabel = "distance from NGRIP [km]", ylabel= "stadial-interstadial d18O difference [o/oo VPDB]",
    framestyle = :box, legend = :topright, xlims=(3500,19000), xticks=(5000:2500:17500),
    formatter = :plain ))
savefig("offsetmedian_distance_cave.pdf")
savefig("offsetmedian_distance_cave.png")


################################################################################
### lat-offset -> caves
plot(df2.mean_off18O[df2.continent.==1],df2.latitude[df2.continent.==1],color = :black,
    seriestype = :scatter, label = "Europe", markersize = 3)
plot!(df2.mean_off18O[df2.continent.==2],df2.latitude[df2.continent.==2],color = :yellow,
    seriestype = :scatter, label = "N America", markersize = 3)
plot!(df2.mean_off18O[df2.continent.==3],df2.latitude[df2.continent.==3], color = :red,
    seriestype = :scatter, label = "S America", markersize = 3)
plot!(df2.mean_off18O[df2.continent.==4],df2.latitude[df2.continent.==4], color = :orange,
    seriestype = :scatter, label = "Africa", markersize = 3)
plot!(df2.mean_off18O[df2.continent.==5],df2.latitude[df2.continent.==5], color = :blue,
    seriestype = :scatter,label = "Asia", markersize = 3)
annotate!(1.9,45, text("b)", color = :black, :left, 12))
display(plot!(df2.mean_off18O[df2.continent.==6],df2.latitude[df2.continent.==6], color = :green,
    seriestype = :scatter, label = "Oceania/Australia", markersize = 3,
    ylabel = "latitude [°]", xlabel= "stadial-interstadial d18O difference  [o/oo VPDB]",
    framestyle = :box, legend = :bottomright))

savefig("offsetmedian_lat_cave.pdf")
savefig("offsetmedian_lat_cave.png")





################################################################################
### GMT codes from here on
################################################################################

# ### test scatter of results with GMT (generic mapping tool)
import GMT
sort!(df2,[:number,:site_id], rev = (true,false))    # order stal after site_id and then after entity_id
lon = df2.longitude[1:size(df2,1)]
lat = df2.latitude[1:size(df2,1)]
off18O = df2.mean_off18O[1:size(df2,1)]
number = df2.number[1:size(df2,1)]
number_print = ones(size(df2,1))*10
number_print = number[number.>9]

### Worldmap
GMT.coast(region=:g, proj=(name=:equidistCylindrical,center=[0 0]), frame=(annot=60, ticks=30, grid=30),
           res=:crude, area=5000, land=:seashell4, water=:white, figsize=14.5) #water=:lightcyan
C = GMT.makecpt(range=(-2.5,2.5,0.5), cmap=:polar)
off18O[off18O.<-2.5] .= -2.49 #for plotting only
GMT.scatter!(lon, lat,zcolor=off18O, markersize=sqrt.(number)/10,#ceil.(number./10).*10 ./100,
    region = [-180 180 -90 90],cmap=C, alpha = 10)#alpha = opacity (25 = 50%))#, show = true
GMT.colorbar!(pos=(paper=true, anchor=(16.3,7), size=(7,0.5), justify=:TC, vertical=true),
    frame=(annot=:auto, ticks=:auto, xlabel="stadial-interstadial @~\144@~@+18@+O difference [\214 VPDB]"),
    fmt = :png, show=true)

### Europe
GMT.coast(region=[-10 40 27 65], proj=(name=:equidistCylindrical,center=[0 0]), frame=(annot=10, ticks=10, grid=10),
           res=:fine, area=5000, land=:lightgrey, water=:white, figsize=14.5) #water=:lightcyan
C = GMT.makecpt(range=(-2.5,2.5,0.5), cmap=:polar)
GMT.scatter!(lon, lat,zcolor=off18O, markersize=sqrt.(number)/5,#ceil.(number./10).*10 ./100,
    region = [-10 40 27 65],cmap=C, alpha = 25)#alpha = opacity (25 = 50%))#, show = true
GMT.colorbar!(pos=(paper=true, anchor=(16.3,9), size=(7,0.5), justify=:TC, vertical=true),
    frame=(annot=:auto, ticks=:auto, xlabel="stadial-interstadial @~\144@~@+18@+O difference [\214 VPDB]"),
    fmt = :png, show=true)#, save = "save.ps")
#GMT.colorbar!(pos=(paper=true, anchor=(16.3,7), size=(7,0.5), justify=:TC, vertical=true),
#    frame=(annot=:auto, ticks=:auto, xlabel="stadial-interstadial @~\144@~@+18@+O [\214 VPDB] offset"),
#    fmt = :png, show=true)
