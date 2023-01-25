function B = Bfunc_free(x,W,w)
B = zeros(size(x));
B(x>=w & x<=W) = 1.0;
f = w./(w-x) - w./x; g = 0.5*(1.0+tanh(f));
B(x<w) = g(x<w);
f = w./(x-W-w) + w./(x-W); g = 0.5*(1.0+tanh(f));
B(x>W & x<W+w) = g(x>W & x<W+w);
end
