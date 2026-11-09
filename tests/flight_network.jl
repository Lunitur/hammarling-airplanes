using Test
using LinearAlgebra
include(joinpath(@__DIR__, "..", "flight_network.jl"))
using .FlightNetwork

@testset "Mreža aerodroma" begin
    mktempdir() do directory
        airports_path = joinpath(directory, "airports.csv")
        flights_path = joinpath(directory, "flights.csv")
        write(airports_path, "IATA_CODE,AIRPORT,LATITUDE,LONGITUDE\nAAA,\"Airport, A\",40,-80\nBBB,Airport B,41,-81\nCCC,Airport C,42,-82\nBAD,Missing coordinates,,\n")
        write(flights_path, "MONTH,ORIGIN_AIRPORT,DESTINATION_AIRPORT,CANCELLED\n1,AAA,BBB,0\n1,AAA,BBB,0\n1,BBB,AAA,0\n1,BBB,CCC,0\n1,AAA,CCC,1\n1,12345,AAA,0\n1,,BBB,0\n1,AAA,AAA,0\n1,AAA,BBB,\n2,AAA,CCC,0\n")
        airports = load_airports(airports_path)
        @test length(airports) == 3
        @test airports["AAA"].name == "Airport, A"
        data = load_routes(flights_path, airports)
        @test (data.selected, data.cancelled, data.invalid, data.accepted) == (9, 1, 4, 4)
        network = build_network(data.routes, airports; n=3)
        @test network.codes == ["BBB", "AAA", "CCC"]
        @test sum(network.counts) == 4
        @test network.W ≈ [0 0.75 0.25; 0.75 0 0; 0.25 0 0]
        @test network.L * ones(3) ≈ zeros(3) atol=1e-14
        @test eigmin(Symmetric(network.L)) >= -1e-14
        @test eigmax(Symmetric(system_matrix(network.L; gamma=0.3))) ≈ -0.3
        B = disturbance_matrix(network.codes, ["AAA"]; sigma=2.0, background=0.1)
        @test B * B' ≈ Diagonal([0.01, 4.01, 0.01])
        @test size(disturbance_matrix(network.codes, ["AAA"]; background=0)) == (3, 1)
        @test_throws ArgumentError load_routes(flights_path, airports; month=13)
        @test_throws ArgumentError system_matrix(network.L; gamma=0)
        @test_throws ArgumentError disturbance_matrix(network.codes, ["UNKNOWN"])
        @test_throws ArgumentError disturbance_matrix(network.codes, ["AAA", "AAA"])
        @test_throws ErrorException ensure_dataset(joinpath(directory, "missing"))
        @test ensure_dataset(directory) == directory
    end
end
