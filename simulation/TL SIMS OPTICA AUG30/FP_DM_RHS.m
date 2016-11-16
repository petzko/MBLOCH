function [Y] = FP_DM_RHS( X )

global dat;

% % % % FP CAVITY

% X(1) = r110,
% X(2) = r330;
% X(3) = r220;
% X(4) = r11p;
% X(5) = r33p;
% X(6) = r22p;
% rRES = rRES;
% X(7) = r130;
% X(8) = r13p;
% X(9) = r13m;
% X(10) = n32p;
% X(11) = n32m;
% X(12) = n12p;
% X(13) = n12m;
% X(14) = U;
% X(15) = V;

W = squeeze(dat.W); INJ = dat.INJ; ULL = dat.ULL; LLL = dat.LLL; RES = dat.RES; DEPOP = dat.DEPOP; G = squeeze(dat.G);

rRES = 1-X(1)-X(2)-X(3);

%%%% POPULATIONS
%time resolved dipole moment ratio.
%%%% POPULATIONS

Y(1) = 1i*dat.O13.*(X(7)-conj(X(7))) +(W(ULL,INJ) + W(ULL,DEPOP)).*X(2) + (W(LLL,INJ)+W(LLL,DEPOP)).*X(3) + (W(RES,INJ)+W(RES,DEPOP)).*rRES- G(INJ).*X(1);
lmInteraction = dat.dipR.*(conj(X(14)).*X(10) + conj(X(15)).*X(11));
Y(2) = 1i*dat.O13.*(conj(X(7)) - X(7)) +1i/2*(lmInteraction-conj(lmInteraction)) +  X(1).*W(INJ,ULL) + X(3).*W(LLL,ULL) + rRES.*W(RES,ULL)- G(ULL).*X(2);
Y(3) = -1i/2*(lmInteraction-conj(lmInteraction)) + X(1).*W(INJ,LLL)+ X(2).*W(ULL,LLL) +rRES.*W(RES,LLL)- G(LLL).*X(3);
Y(4) = 1i*dat.O13.*(X(8)-conj(X(9))) + (W(ULL,INJ)+W(ULL,DEPOP)).*X(5)+ (W(LLL,INJ)+W(LLL,DEPOP)).*X(6) - (G(INJ)+dat.diffusion).*X(4);
%%% r33+
Y(5) = 1i*dat.O13.*(conj(X(9))-X(8))+1i/2*dat.dipR.*(conj(X(15)).*(X(10)) -(X(14)).*conj(X(11))) +  W(INJ,ULL).*X(4) + W(LLL,ULL).*X(6) - (G(ULL)+dat.diffusion).*X(5);
%%% r22+
Y(6) = -1i/2*dat.dipR.*(conj(X(15)).*(X(10)) -(X(14)).*conj(X(11))) + W(INJ,LLL).*X(4) + W(ULL,LLL).*X(5) - (G(LLL)+dat.diffusion).*X(6);

%%%% COHERENCES
Y(7) = dat.dE13.*X(7) + 1i*dat.O13.*(X(1)-X(2)) +1i/2*dat.dipR.*(conj(X(14)).*X(12) + conj(X(15)).*X(13));
%%% r13+
Y(8) = (dat.dE13-dat.diffusion).*X(8) + 1i*dat.O13.*(X(4)-X(5)) +1i/2*dat.dipR.*(conj(X(15)).*X(12));
%%% r13-
Y(9) = (dat.dE13-dat.diffusion).*X(9) + 1i*dat.O13.*conj(X(4)-X(5)) + 1i/2*dat.dipR.*(conj(X(14)).*X(13));

Y(10) = dat.dE32.*X(10) + 1i/2*dat.dipR.*(X(14).*(X(2)-X(3)) + X(15).*(X(5)-X(6))) - 1i*dat.O13.*X(12);
Y(11) = dat.dE32.*X(11) + 1i/2*dat.dipR.*(X(15).*(X(2)-X(3)) + X(14).*conj(X(5)-X(6))) - 1i*dat.O13.*X(13);
Y(12) = dat.dE12.*X(12) +1i/2*dat.dipR.*(X(14).*X(7) + X(15).*X(8)) - 1i*dat.O13.*X(10);
Y(13) = dat.dE12.*X(13) + 1i/2*dat.dipR.*(X(15).*X(7) +X(14).*X(9)) - 1i*dat.O13.*X(11);
%FIELDS
factor = dat.factor*dat.dipR;
Y(14)  = factor.*X(10)+dat.losses*X(14);
Y(15)  = factor.*X(11)+dat.losses*X(15);


end