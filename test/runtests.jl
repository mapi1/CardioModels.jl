using Main.CardioModels
using Test


N = 100

@testset "Baselli" begin
    @test_throws DomainError getModel(5)
    model = getModel(1)
    result = predict(model, N)
    @test length(result) == 3
    @test length(result[1]) == N
end

@testset "DeBoer" begin
    model = DeBoerModel()
    result = predict!(model, N)
    @test length(result) == 5
    @test length(result[1]) == N
end

@testset "Karemaker" begin
    model = KaremakerModel()
    result = predict(model, N)
    @test length(result) == 6
    @test length(result[1]) == N
end