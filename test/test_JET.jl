# `JET.test_package` hangs indefinitely on Julia 1.12: JET's LoweredCodeUtils-based
# definition scan hits an `UndefRefError` on a keyword-argument sorter (emitting
# "skipping callee ... due to UndefRefError()" from LoweredCodeUtils/signatures.jl)
# and then stalls. The package's runtime code itself analyzes/tests cleanly, so this
# is a JET/LoweredCodeUtils-vs-Julia-1.12 limitation. Run the static analysis on
# 1.10/1.11 and skip it on 1.12+ until the upstream issue is resolved.
const _SKIP_JET = VERSION == v"1.12"

if _SKIP_JET
    @info "Skipping JET.test_package on Julia $(VERSION): JET/LoweredCodeUtils definition scan hangs on 1.12+."
else
    @testset "JET Testing" begin
        rep = JET.test_package(
            AstroForceModels;
            toplevel_logger=nothing,
            target_modules=(@__MODULE__,),
            analyze_from_definitions=false, #TODO: REMOVE THIS LATER, SOMETHING IN ALBEDO CAUSE THIS TO TAKE HOURS, BUT ALBEDO DOES EVENUTALLY PASS
        )
    end
end
