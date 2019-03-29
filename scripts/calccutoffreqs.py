from oct2py import octave

def cutoffreqs(M, q):
  
  octave.eval("clear") 
  octave.eval("addpath('../src');")
  octave.eval("setconstants") 

  octave.push("m", M);
  octave.push("q", q);

  octave.eval("eta = q/(1+q)^2;") 

  # calculate the mass and spin of the final BH 
  octave.eval("[mf, af, A, Omega] = finalmassandspin_eobnrv2(m, eta);")

  # QNM freq (l = 2, m = 2, n = 0 mode) - Berti's formula 
  octave.eval("f_QNM = real(Omega(1))/(2*pi*m*LAL_MTSUN_SI);")
  octave.eval("tau = -1/imag(Omega(1))*m*LAL_MTSUN_SI;")
  octave.eval("Q = tau*f_QNM*sqrt(2)*pi;")
  octave.eval("sigma = f_QNM/Q;")

  # ISCO freq - Schwarzschild 
  octave.eval("[flso_Schw, flRing_Schw] = calcflso(m);")

  # ISCO freq - Kerr. Also QNM freq using old formulas -- for cross checking 
  octave.eval("[flso_Kerr, flRing_Kerr, fQNM_old, Q_old] = CalcFlsoKerrBH(mf, af);")
  
  flso_Kerr = octave.pull("flso_Kerr")
  flRing_Kerr = octave.pull("flRing_Kerr")
  
  return (flso_Kerr, flRing_Kerr)
