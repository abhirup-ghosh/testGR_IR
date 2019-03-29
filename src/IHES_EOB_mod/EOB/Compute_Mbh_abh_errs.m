function err = Compute_Mbh_abh_errs(nqnm,lmaxqnm, EOBopt,Mbh,abh,Omega,t, dt, wav,nu,H0,pph0)
%Compute_Mbh_abh_errs Computes the norm of the combined errors in Mbh and abh from the given
%values compared to those obtained from the waveform using energy and angular momentum
%balance

sigma = EOBQNMKerrFit_mod(nqnm(1),lmaxqnm, EOBopt.L,EOBopt.M,abh, EOBopt.QNMDataDir, EOBopt.iQNMfac);
if nqnm(2)>0
    sigma = [sigma conj(sigma(1:nqnm(2)))];
end
sigma = sigma/Mbh; % rescale units !
[psilm,dpsilm,omglm] = EOBhlmQNMMatch(Omega,t, sigma, nqnm, dt, wav );

wav        = EOBExtractAPO(psilm,t);  
wav.t      = t;
wav.psilm  = psilm;

Erad = 0;
Jrad = 0;

LM2K = EOBopt.LM2K;

% Multiply the 2,1, 3,3, and 3,1 psilms and dlms by a1^0.5

wav.psilm(LM2K(2,1),:) = wav.psilm(LM2K(2,1),:) .* EOBopt.a1^0.5;
wav.psilm(LM2K(3,3),:) = wav.psilm(LM2K(3,3),:) .* EOBopt.a1^0.5;
wav.psilm(LM2K(3,1),:) = wav.psilm(LM2K(3,1),:) .* EOBopt.a1^0.5;

wav.dlm(LM2K(2,1),:) = wav.dlm(LM2K(2,1),:) .* EOBopt.a1^0.5;
wav.dlm(LM2K(3,3),:) = wav.dlm(LM2K(3,3),:) .* EOBopt.a1^0.5;
wav.dlm(LM2K(3,1),:) = wav.dlm(LM2K(3,1),:) .* EOBopt.a1^0.5;

% Multiply the 3,2, 4,4, and 4,2 psilms and dlms by a2^0.5

wav.psilm(LM2K(3,2),:) = wav.psilm(LM2K(3,2),:) .* EOBopt.a2^0.5;
wav.psilm(LM2K(4,4),:) = wav.psilm(LM2K(4,4),:) .* EOBopt.a2^0.5;
wav.psilm(LM2K(4,2),:) = wav.psilm(LM2K(4,2),:) .* EOBopt.a2^0.5;

wav.dlm(LM2K(3,2),:) = wav.dlm(LM2K(3,2),:) .* EOBopt.a2^0.5;
wav.dlm(LM2K(4,4),:) = wav.dlm(LM2K(4,4),:) .* EOBopt.a2^0.5;
wav.dlm(LM2K(4,2),:) = wav.dlm(LM2K(4,2),:) .* EOBopt.a2^0.5;

for ll = 2:7 %(formerly we took the sum to ll = 8; we now restrict to ll = 7 to match the modes we can fit QNMs to)
    for mm = 1:ll
        Erad = Erad + 2*EOBopt.fac*(ll - 1)*ll*(ll + 1)*(ll + 2)*nu^2*trapz(wav.t,abs(wav.dlm(LM2K(ll,mm),:)).^2)/(16*pi);
        Jrad = Jrad + 2*EOBopt.fac*(ll - 1)*ll*(ll + 1)*(ll + 2)*mm*nu^2*trapz(wav.t,imag(wav.psilm(LM2K(ll,mm),:).*conj(wav.dlm(LM2K(ll,mm),:))))/(16*pi);
    end
end

Mbhrad = nu*H0 - Erad;
abhrad = (nu*pph0 - Jrad)/Mbhrad^2;

err = sqrt((Mbh/Mbhrad-1)^2+(abh/abhrad-1)^2);

end

