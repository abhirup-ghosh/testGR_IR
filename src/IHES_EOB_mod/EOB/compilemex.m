% matlab script to compile mex-files

% linux
mex Lagrange1DIntMex.c -v CFLAGS="\$CFLAGS -Wall" -output LagInt1d
mex ModTailMex.c -v CFLAGS="\$CFLAGS -Wall" -output EOBModTailm
%mex GammaComplexMex.c -v CFLAGS="\$CFLAGS -Wall" -o GammaComplexm