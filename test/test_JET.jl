@testset "JET Testing" begin
    rep = JET.test_package(
        AstroForceModels;
        toplevel_logger=nothing,
        target_modules=(@__MODULE__,),
        analyze_from_definitions=false, #TODO: REMOVE THIS LATER, SOMETHING IN ALBEDO CAUSE THIS TO TAKE HOURS, BUT ALBEDO DOES EVENUTALLY PASS
    )
end
