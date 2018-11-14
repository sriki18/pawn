# PAWN
This is a MATLAB implementation of the global sensitivity analysis algorithm PAWN.
- `PAWN.m` contains the MATLAB implementation of the PAWN algorithm.
- `ishigami_homma.m` recreates Fig. 4 of the paper [1].

## Example

```matlab
% Ishigami-Homma function
f = @(x,p) sin(x(1)) + p.a * sin(x(2)^2) + p.b * x(3)^4 * sin(x(1));
p.a = 2; p.b = 1;
ih = @(x) f(x, p);

% Bounds
lb = ones(1,3) * -pi; ub = ones(1,3) * pi;

% Parameters from Figure 4's caption
n = 15;
Nu = 100;
Nc = 100;

% Other parameters
npts = 100; seed = 4;

[KS,xvals,y_u, y_c, par_u, par_c, ft] = PAWN(f, p, lb, ub, Nu, n, Nc, npts, seed);
```

## References
[1] Pianosi, F., Wagener, T., 2015. A simple and efficient method for global sensitivity analysis based on cumulative distribution functions. Environ. Model. Softw. 67, 1–11. www.doi.org/10.1016/j.envsoft.2015.01.004