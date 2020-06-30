function [poly_bezier,dpoly_bezier] = poly_eval_bezier(s,M,a)

poly_bezier  = 0;
dpoly_bezier = 0;

for k = 0 : M
    val = factorial(M)/( factorial(k)*factorial(M-k) );
    term = a(k+1) * val * s^k * (1-s)^(M-k);
    poly_bezier = poly_bezier + term;
end

if M == 0
    dpoly_bezier = 0;
else
    for k = 0 : M-1
        val = factorial(M)/factorial(k)/factorial(M-k-1);
        term  = (a(k+2) - a(k+1)) * val * s^k * (1-s)^(M-k-1);
        dpoly_bezier = dpoly_bezier + term;
    end
end