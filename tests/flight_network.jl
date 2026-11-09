using Test
using Dates
using LinearAlgebra
using Parquet2
using Tables
include(joinpath(@__DIR__, "..", "flight_network.jl"))
using .FlightNetwork

@testset "Europska mreža aerodroma (OPDI)" begin
    mktempdir() do directory
        airports_path = joinpath(directory, "airports.csv")
        flights_path = joinpath(directory, "flight_list_202501.parquet")
        write(airports_path, """
        ident,icao_code,name,latitude_deg,longitude_deg,continent,type
        LOCAL,EAAA,"Aerodrom, čvor A",45,16,EU,large_airport
        EBBB,EBBB,Aerodrom B,48,11,EU,large_airport
        ECCC,ECCC,Aerodrom C,51,0,EU,medium_airport
        EBAD,EBAD,Bez koordinata,,,EU,small_airport
        ENAN,ENAN,Nevaljane koordinate,NaN,10,EU,small_airport
        EOUT,EOUT,Izvan raspona,91,10,EU,small_airport
        ECLO,ECLO,Zatvoren,45,16,EU,closed
        KJFK,KJFK,New York,40,-73,NA,large_airport
        XXXX,,Bez ICAO oznake,45,16,EU,small_airport
        """)
        flights = (
            dof=Union{Missing, DateTime}[fill(DateTime(2025, 1, 15), 9);
                DateTime(2025, 2, 1); DateTime(2024, 1, 15); missing],
            adep=Union{Missing, String}["EAAA", "EAAA", "EBBB", "EBBB", "KJFK",
                missing, "EAAA", "EAAA", "", "EAAA", "EAAA", "EAAA"],
            ades=Union{Missing, String}["EBBB", "EBBB", "EAAA", "ECCC", "EAAA",
                "EBBB", "EAAA", missing, "EBBB", "ECCC", "ECCC", "ECCC"],
            adep_p=fill("EAAA", 12),
            ades_p=fill("ECCC", 12),
        )
        parts = [map(col -> col[1:6], flights), map(col -> col[7:end], flights)]
        Parquet2.writefile(flights_path, Tables.partitioner(parts))
        @test length(Parquet2.Dataset(flights_path)) == 2
        airports = load_airports(airports_path)
        @test sort(collect(keys(airports))) == ["EAAA", "EBBB", "ECCC"]
        @test airports["EAAA"].name == "Aerodrom, čvor A"
        data = load_routes(flights_path, airports)
        @test (data.selected, data.outside_period, data.invalid, data.accepted) == (9, 3, 5, 4)
        @test data.selected == data.invalid + data.accepted
        @test sum(values(data.routes)) == data.accepted
        @test data.routes == Dict(("EAAA", "EBBB") => 2, ("EBBB", "EAAA") => 1, ("EBBB", "ECCC") => 1)
        @test load_routes(flights_path, airports; year=2024).accepted == 1
        @test load_routes(flights_path, airports; month=2).accepted == 1
        network = build_network(data.routes, airports; n=3)
        @test network.codes == ["EBBB", "EAAA", "ECCC"]
        @test sum(network.counts) == 4
        @test network.W ≈ [0 0.75 0.25; 0.75 0 0; 0.25 0 0]
        @test network.L * ones(3) ≈ zeros(3) atol=1e-14
        @test eigmin(Symmetric(network.L)) >= -1e-14
        @test eigmax(Symmetric(system_matrix(network.L; gamma=0.3))) ≈ -0.3
        B = disturbance_matrix(network.codes, ["EAAA"]; sigma=2.0, background=0.1)
        @test B * B' ≈ Diagonal([0.01, 4.01, 0.01])
        @test size(disturbance_matrix(network.codes, ["EAAA"]; background=0)) == (3, 1)
        @test_throws ArgumentError load_routes(flights_path, airports; month=13)
        @test_throws ArgumentError load_routes(flights_path, airports; year=2021)
        @test_throws ErrorException load_routes(flights_path, airports; month=3)
        @test_throws ArgumentError system_matrix(network.L; gamma=0)
        @test_throws ArgumentError disturbance_matrix(network.codes, ["UNKNOWN"])
        @test_throws ArgumentError disturbance_matrix(network.codes, ["EAAA", "EAAA"])
        @test_throws ErrorException ensure_dataset(joinpath(directory, "missing"))
        @test_throws ErrorException ensure_dataset(directory; month=2)
        @test_throws ArgumentError ensure_dataset(directory; month=0)
        paths = (flights=flights_path, airports=airports_path)
        @test ensure_dataset(directory) == paths
        @test ensure_dataset(directory; download=true) == paths

        bad_flights = joinpath(directory, "bad.parquet")
        Parquet2.writefile(bad_flights, (adep=["EAAA"], ades=["EBBB"]))
        @test_throws ErrorException load_routes(bad_flights, airports)
        write(airports_path, "IATA_CODE,AIRPORT,LATITUDE,LONGITUDE\nAAA,Legacy,45,16\n")
        @test_throws ArgumentError load_airports(airports_path)
    end
end
