

# Testing the polynomial evaluation macro:
p = randn(6)
x = 0.33
@test @evalpoly(x, p[1], p[2]) ≈ Psychro.@polyeval(x, p, 2) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3]) ≈ Psychro.@polyeval(x, p, 3) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4]) ≈ Psychro.@polyeval(x, p, 4) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4], p[5]) ≈ Psychro.@polyeval(x, p, 5) atol=1e-13
@test @evalpoly(x, p[1], p[2], p[3], p[4], p[5], p[6]) ≈ Psychro.@polyeval(x, p, 6) atol=1e-13
