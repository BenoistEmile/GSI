#include "Model.h"
#include <julia.h>
#include "fmt/format.h"
#include "juliasolver.h"

JuliaSolver::JuliaSolver()
{
  jl_init();
  jl_eval_string("include(\"./ juliasolver.jl \")");
  jl_eval_string("import .juliasolver as jl");
  jl_eval_string("jl.init()");
}

JuliaSolver::~JuliaSolver()
{
  jl_atexit_hook(0);
}

void
JuliaSolver::addQVars(std::size_t size)
{
  m_num_Q = size;
  jl_eval_string(fmt::format("@variable(model, P[1:{0:.0f}])", size));
  jl_eval_string(fmt::format("@variable(model, Q[1:{0:.0f}])", size));
}

void
JuliaSolver::addXVars(std::size_t size)
{
  m_num_X = size;
}

void (*func_jl)() = jl_unbox_voidpointer(jl_eval_string("@cfunction()"))