ENV["GENX_PRECOMPILE"] = "false"
using GenX
GenX.run_genx_case!(dirname(@__FILE__))
