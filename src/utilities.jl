"""
    @polyeval(z, p, N)

Evaluate the polynomial ``\\sum_k c[k] z^{k-1}`` for the coefficients `c[1]`, `c[2]`, ..., `c[N]`;
that is, the coefficients are given in ascending order by power of `z`.  This macro expands
to efficient inline code that uses either Horner's method.

```jldoctest

julia> p = randn(5)
5-element Array{Float64,1}:
  0.106455 
  0.0351716
 -0.204764 
 -0.87642  
  0.260407 

julia> @polyeval(0.5, p, 5)
-0.020427265391245605

```
"""
macro polyeval(x, p, N)
    if N==1
        return :($(esc(p))[1])
    elseif N==2
        return :($(esc(p))[1] + $(esc(x))*$(esc(p))[2])
    else
        nm1 = N-1
        ex = :($(esc(p))[$nm1] + $(esc(x))*$(esc(p))[$N])
        for i in N-2:-1:1
            ex = :($(esc(p))[$i] + $(esc(x))*$ex)
        end
    end
    return ex
end

"""
    ConvergenceError([msg, val, niter, err])

Exception that is to be thrown when an iterative procedure fails
to converge.

"""
mutable struct ConvergenceError <: Exception
    msg::String
    val::Float64
    niter::Int
    err::Float64
    ConvergenceError() = new("Calculations failed to converge", NaN, -1, NaN)
    ConvergenceError(msg) = new(msg, NaN, -1, NaN)
    ConvergenceError(msg, val) = new(msg, val, -1, NaN)
    ConvergenceError(msg, val, niter) = new(msg, val, niter, NaN)
    ConvergenceError(msg, val, niter, err) = new(msg, val, niter, err)
end
