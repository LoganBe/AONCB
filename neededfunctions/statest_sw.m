function [mu,sig2] = statest_sw(theta,rhoe,rhoi,rhoei)

if nargin == 1; rhoe = 0; rhoi = 0; rhoei = 0; end

re = theta.r.ee; Re = re*1e-3;
ri = theta.r.ei; Ri = ri*1e-3;
Ke = theta.K.ee;
Ki = theta.K.ei;
Ve = theta.ereversal;
Vi = theta.ireversal;
tauleak = theta.tauleak;
W = theta.Wf;
if ~iscell(W)
    We = W(1,1); Wi = W(1,2);
else
    We = W{1,1}; Wi = W{1,2};
end
if ~isfield(theta,'I'); I = 0; else; I = theta.I; end
%if length(re) > 1; [Re,Ri] = meshgrid(re,ri); else; Re = re; Ri = ri; end
%if length(Ke) > 1; [Ke,Ki] = meshgrid(Ke,Ki); e
if ~isfield(theta,'ir'); irflag = false; else; irflag = theta.ir; end %Default no intrinisc resistance

if irflag
    if ~isfield(theta,'meanv'); error("Need Desired Mean V"); end
    mv = theta.meanv;
    gv = theta.irf;
    gv = 1+gv(mv);
else
    gv = 1;
end

mu = (Ke*Re*We*Ve + Ki*Ri*Wi*Vi + I)./(gv + Ke*Re*We + Ki*Ri*Wi);

sig2 =  (1./(2*tauleak*(gv + Ke.*Re.*We + Ki.*Ri.*Wi))).* ...
       (Ke.*Re.*We.^2.*(1+rhoe.*(Ke - 1)).*(Ve - mu).^2 + ...
        Ki.*Ri.*Wi.^2.*(1+rhoi.*(Ki - 1)).*(Vi - mu).^2 - ...
        2.*rhoei.*sqrt(Re.*Ri).*Ke.*We.*Ki.*Wi.*(Ve - mu).*(mu - Vi));