a = 0;
b = 1;
epsr = 1e-05;
nlimit = int64(0);
[result, npts, relerr, ifail] = d01ah(a, b, epsr, 'sin', nlimit)