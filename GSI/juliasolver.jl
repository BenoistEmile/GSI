module juliasolver

using JuMP
# using HiGHS
using CPLEX
import MultiObjectiveAlgorithms as MOA

function init()
    global model = JuMP.Model(() -> MOA.Optimizer(CPLEX.Optimizer))
    set_attribute(model, "CPX_PARAM_EPINT", 1e-8)
    set_attribute(model, MOA.Algorithm(), MOA.DominguezRios())
    set_attribute(model, MOI.TimeLimitSec(), 20)
end

function test()
    @variable(model, P1, Bin)
    @variable(model, P2, Bin)
end

end