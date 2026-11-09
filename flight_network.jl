module FlightNetwork

using CSV
using Dates
using Downloads
using LinearAlgebra
using Parquet2
using Tables

export ensure_dataset, load_airports, load_routes, build_network, system_matrix,
       disturbance_matrix, airport_map, covariance_heatmap, recovery_plot, result_table

const AIRPORTS_URL = "https://raw.githubusercontent.com/davidmegginson/ourairports-data/7dd521a50a421171db1cb69f87657a60b48561c7/airports.csv"

function flight_period(year, month)
    year >= 2022 || throw(ArgumentError("OPDI podaci počinju 2022. godine."))
    month in 1:12 || throw(ArgumentError("Mjesec mora biti između 1 i 12."))
    Date(year, month, 1)
end

function ensure_dataset(directory; year=2025, month=1, download=false)
    period = Dates.format(flight_period(year, month), "yyyymm")
    filename = "flight_list_$period.parquet"
    paths = (flights=joinpath(directory, filename), airports=joinpath(directory, "airports.csv"))
    all(isfile, paths) && return paths
    download || error("Nedostaju $filename ili OurAirports airports.csv u $directory. " *
        "Preuzmite podatke s https://www.opdi.aero/flight-list-data i $AIRPORTS_URL " *
        "ili postavite download_data = true u bilježnici.")
    mkpath(directory)
    urls = ("https://www.eurocontrol.int/performance/data/download/OPDI/v002/flight_list/$filename", AIRPORTS_URL)
    for (path, url) in zip(paths, urls)
        isfile(path) && continue
        try
            mktempdir(directory) do temporary
                target = joinpath(temporary, basename(path))
                Downloads.download(url, target)
                if endswith(path, ".parquet")
                    flight_columns(Parquet2.Dataset(target))
                else
                    load_airports(target)
                end
                mv(target, path; force=true)
            end
        catch exception
            error("Preuzimanje nije uspjelo: $(sprint(showerror, exception)). " *
                "Preuzmite $url ručno u $path.")
        end
    end
    paths
end

function load_airports(path)
    columns = [:icao_code, :name, :latitude_deg, :longitude_deg, :continent, :type]
    rows = CSV.File(path; select=columns,
        types=Dict(:icao_code => String, :name => String, :continent => String, :type => String,
                   :latitude_deg => Float64, :longitude_deg => Float64), strict=true)
    all(name -> name in propertynames(rows), columns) || error("Neispravni OurAirports stupci u airports.csv.")
    airports = Dict{String, NamedTuple{(:code, :name, :lat, :lon), Tuple{String, String, Float64, Float64}}}()
    for row in rows
        isequal(row.continent, "EU") && !ismissing(row.type) && row.type != "closed" || continue
        any(ismissing, (row.icao_code, row.name, row.latitude_deg, row.longitude_deg)) && continue
        isfinite(row.latitude_deg) && isfinite(row.longitude_deg) || continue
        abs(row.latitude_deg) <= 90 && abs(row.longitude_deg) <= 180 || continue
        code = String(row.icao_code)
        occursin(r"^[A-Z]{4}$", code) || continue
        haskey(airports, code) && error("Ponovljena ICAO oznaka u airports.csv: $code.")
        airports[code] = (code=code, name=String(row.name), lat=row.latitude_deg, lon=row.longitude_deg)
    end
    isempty(airports) && error("Nema europskih aerodroma s ICAO oznakama i valjanim koordinatama.")
    airports
end

function flight_columns(dataset)
    columns = (:dof, :adep, :ades)
    all(name -> name in Tables.schema(dataset).names, columns) ||
        error("OPDI Parquet mora sadržavati stupce dof, adep i ades.")
    columns
end

function load_routes(path, airports; year=2025, month=1)
    start = flight_period(year, month)
    stop = start + Month(1)
    dataset = Parquet2.Dataset(path)
    columns = flight_columns(dataset)
    routes = Dict{Tuple{String, String}, Int}()
    selected, outside_period, invalid, accepted = 0, 0, 0, 0
    for group in dataset
        dates, origins, destinations = map(name -> Tables.getcolumn(group, name), columns)
        for (date, origin, destination) in zip(dates, origins, destinations)
            if ismissing(date) || !(start <= Date(date) < stop)
                outside_period += 1
                continue
            end
            selected += 1
            if ismissing(origin) || ismissing(destination) ||
               !haskey(airports, origin) || !haskey(airports, destination) || origin == destination
                invalid += 1
                continue
            end
            key = (String(origin), String(destination))
            routes[key] = get(routes, key, 0) + 1
            accepted += 1
        end
    end
    isempty(routes) && error("Nema upotrebljivih letova za $year-$month između europskih aerodroma.")
    (routes=routes, selected=selected, outside_period=outside_period, invalid=invalid, accepted=accepted)
end

function build_network(routes, airports; n=50)
    n >= 2 || throw(ArgumentError("Odaberite barem dva aerodroma."))
    traffic = Dict{String, Int}()
    for ((origin, destination), count) in routes
        count > 0 || throw(ArgumentError("Broj letova mora biti pozitivan."))
        for code in (origin, destination)
            traffic[code] = get(traffic, code, 0) + count
        end
    end
    codes = sort(collect(keys(traffic)); by=code -> (-traffic[code], code))
    length(codes) >= 2 || error("Mreža treba barem dva aerodroma.")
    codes = codes[1:min(n, length(codes))]
    indices = Dict(code => i for (i, code) in enumerate(codes))
    counts = zeros(length(codes), length(codes))
    for ((origin, destination), count) in routes
        haskey(indices, origin) && haskey(indices, destination) || continue
        counts[indices[origin], indices[destination]] += count
    end
    weights = (counts + counts') / 2
    scale = maximum(vec(sum(weights; dims=2)))
    scale > 0 || error("Odabrani aerodromi nemaju međusobnih veza. Povećajte n.")
    W = weights / scale
    L = Matrix(Diagonal(vec(sum(W; dims=2))) - W)
    (codes=codes, airports=[airports[code] for code in codes], counts=counts,
     W=W, L=L, traffic=[traffic[code] for code in codes], weight_scale=scale)
end

function system_matrix(L; alpha=1.0, gamma=0.2)
    isfinite(alpha) && alpha > 0 && isfinite(gamma) && gamma > 0 ||
        throw(ArgumentError("alpha i gamma moraju biti konačni i pozitivni."))
    -alpha * L - gamma * I
end

function disturbance_matrix(codes, sources; sigma=1.0, background=0.01)
    isfinite(sigma) && sigma > 0 && isfinite(background) && background >= 0 ||
        throw(ArgumentError("sigma mora biti pozitivna, a background nenegativan."))
    isempty(sources) && throw(ArgumentError("Odaberite barem jedan izvor."))
    length(unique(sources)) == length(sources) || throw(ArgumentError("Izvori se ponavljaju."))
    indices = [findfirst(==(source), codes) for source in sources]
    any(isnothing, indices) && throw(ArgumentError("Izvor nije u odabranoj mreži."))
    E = Matrix{Float64}(I, length(codes), length(codes))
    localized = sigma * E[:, Int.(indices)]
    iszero(background) ? localized : hcat(localized, background * E)
end

escape_html(value) = replace(string(value), '&' => "&amp;", '<' => "&lt;", '>' => "&gt;",
    '"' => "&quot;", '\'' => "&#39;")
color(t) = "rgb($(round(Int, 235 - 205clamp(t, 0, 1))),$(round(Int, 245 - 125clamp(t, 0, 1))),$(round(Int, 255 - 90clamp(t, 0, 1))))"

function result_table(rows)
    isempty(rows) && return Base.HTML("<p>Nema rezultata.</p>")
    names = propertynames(first(rows))
    header = join(("<th>$(escape_html(name))</th>" for name in names))
    body = join(("<tr>" * join(("<td>$(escape_html(getproperty(row, name)))</td>" for name in names)) * "</tr>" for row in rows))
    Base.HTML("<table><thead><tr>$header</tr></thead><tbody>$body</tbody></table>")
end

function airport_map(network, variance)
    nodes, W = network.airports, network.W
    lon, lat = getproperty.(nodes, :lon), getproperty.(nodes, :lat)
    project(values, offset, span) = offset .+ span .* (values .- minimum(values)) ./ max(maximum(values) - minimum(values), 1)
    x, y = project(lon, 40, 820), 450 .- project(lat, 0, 390)
    maximum_variance = max(maximum(variance), eps())
    shapes = String[]
    for j in 2:length(nodes), i in 1:j-1
        W[i, j] > 0 || continue
        push!(shapes, "<line x1='$(x[i])' y1='$(y[i])' x2='$(x[j])' y2='$(y[j])' stroke='#94a3b8' stroke-opacity='0.25' stroke-width='$(0.3 + 3W[i, j])'/>")
    end
    for i in eachindex(nodes)
        tooltip = escape_html("$(nodes[i].code): $(nodes[i].name), Xᵢᵢ=$(round(variance[i]; sigdigits=4))")
        push!(shapes, "<circle cx='$(x[i])' cy='$(y[i])' r='$(4 + 10sqrt(max(variance[i], 0) / maximum_variance))' fill='$(color(variance[i] / maximum_variance))' stroke='#334155'><title>$tooltip</title></circle>")
        push!(shapes, "<text x='$(x[i]+7)' y='$(y[i]-7)' font-size='10'>$(escape_html(nodes[i].code))</text>")
    end
    Base.HTML("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 920 500' style='max-width:100%;background:#fff;color:#111;font-family:sans-serif'><text x='20' y='25'>Mreža aerodroma: tamnije / veće = veća modelirana varijanca</text>$(join(shapes))<text x='20' y='485' font-size='12'>Geografski raspored (dužina / širina), bez kartografske podloge. Detalji: prijeđite mišem preko čvora.</text></svg>")
end

function covariance_heatmap(codes, X; limit=20)
    n = min(limit, length(codes))
    scale = max(maximum(abs, X[1:n, 1:n]), eps())
    shapes = String[]
    for j in 1:n, i in 1:n
        value = X[i, j]
        push!(shapes, "<rect x='$(90+22(j-1))' y='$(65+22(i-1))' width='22' height='22' fill='$(color(abs(value)/scale))'><title>$(escape_html(codes[i])) / $(escape_html(codes[j])): $value</title></rect>")
    end
    for i in 1:n
        label = escape_html(codes[i])
        push!(shapes, "<text x='80' y='$(80+22(i-1))' text-anchor='end' font-size='10'>$label</text><text transform='translate($(105+22(i-1)),55) rotate(-55)' font-size='10'>$label</text>")
    end
    Base.HTML("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 $(120+22n) $(90+22n)' style='max-width:650px;width:100%;background:#fff;font-family:sans-serif'><text x='15' y='18'>|Xᵢⱼ| za prvih $n aerodroma; tamnije = veće</text>$(join(shapes))</svg>")
end

function recovery_plot(rows)
    gammas, totals = getproperty.(rows, :gamma), getproperty.(rows, :total_variance)
    low, high = extrema(log10.(gammas))
    x = 70 .+ 610 .* (log10.(gammas) .- low) ./ max(high-low, eps())
    y = 320 .- 260 .* totals ./ max(maximum(totals), eps())
    points = join(("$(x[i]),$(y[i])" for i in eachindex(rows)), " ")
    labels = join(("<circle cx='$(x[i])' cy='$(y[i])' r='4' fill='#155e75'><title>γ=$(gammas[i]), tr(X)=$(totals[i])</title></circle><text x='$(x[i])' y='340' text-anchor='middle' font-size='11'>$(gammas[i])</text>" for i in eachindex(rows)))
    Base.HTML("<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 740 390' style='max-width:100%;background:#fff;font-family:sans-serif'><text x='30' y='25'>Ukupna stacionarna varijanca tr(X) i brzina oporavka γ</text><path d='M70 50 V320 H700' fill='none' stroke='#475569'/><text x='60' y='324' text-anchor='end'>0</text><text x='60' y='64' text-anchor='end' font-size='11'>$(round(maximum(totals); sigdigits=3))</text><polyline points='$points' fill='none' stroke='#155e75' stroke-width='2'/>$labels<text x='300' y='375' font-size='13'>γ (logaritamska os)</text></svg>")
end

end
