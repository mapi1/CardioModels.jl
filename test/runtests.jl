using CardioModels
using Test


N = 100

@testset "Baselli" begin
    @test_throws DomainError getModel(5)
    model = getModel(1)
    result = predict(model, N)
    @test length(result) == 3
    @test length(result[1]) == N

    model_est = fitBaselli(result...)
    @test typeof(model_est) == BaselliModel

    post = postprocess(model_est, result[1], result[2])
    @test length(post) == 8
end

@testset "DeBoer" begin
    model = DeBoerModel()
    result = predict(model, N)
    @test length(result) == 5
    @test length(result[1]) == N

    result = predict!(model, N)
    @test model.hasState
    
    result = predict!(model, 1)
    @test length(result[1]) == 1
    @test model.S[end] == result[1][1]

end

@testset "Karemaker" begin
    model = KaremakerModel()
    result = predict(model, N)
    @test length(result) == 6
    @test length(result[1]) == N
end