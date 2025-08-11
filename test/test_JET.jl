@testset "JET Testing" begin
    rep = JET.test_package(
        AstroForceModels; toplevel_logger=nothing, target_modules=(@__MODULE__,)
    )
end
