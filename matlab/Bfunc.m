function B = Bfunc(x,W,w)
xa = abs(x);
B = zeros(size(x));
B(xa<=W) = 1.0;
f = w./(xa-W-w) + w./(xa-W);
g = 0.5*(1.0+tanh(f));
B(xa>W & xa<W+w) = g(xa>W & xa<W+w);
end
