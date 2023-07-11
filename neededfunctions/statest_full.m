function[meanV, varV, outinfo] = statest_full(theta)

%STATCORR_BETA: Theoertical Mean and Variance of V
%[meanV, varV] = statcorr_beta(theta) Computes the mean and variance of an
%                AONCB neuron. If input correlations are present, assumes a
%                Beta distribution for the deFinetti Measure
%
% Input:
%       theta -- structure with model paratmers given by:
%                theta.K.ee = #E->E Synapses. Scalar
%                theta.r.ee = Excitatory Firing rate (Hz). Scalar
%                theta.ereversal = Excitatory reversal potential (mV). Scalar 
%                theta.K.ei = #I->E Synapses. Scalar
%                theta.r.ei = Inhibitory Firing rate (Hz). Scalar
%                theta.ireversal = Inhibitory reversal potential (mV). Scalar
%                theta.Wf = Synaptic weights for [E->E,I->E; E->I, I->I]. [2x2] mat
%                theta.tauleak = Leak time constant (ms). Scalar
%                theta.Ie = External Current (A). Scalar (If not present, assumes = 0)
%                theta.corrinfo.ae1 = alpha (excitatory) paramter for Beta(). Scalar
%                theta.corrinfo.ae2 = beta (excitatory) paramter for Beta(). Scalar
%                theta.corrinfo.ai1 = alpha (inhibitory) paramter for Beta(). Scalar
%                theta.corrinfo.ai2 = beta (inhibitory) paramter for Beta(). Scalar
%                theta.corrinfo.corridx = arho (Determins if there is cross corr. Assumes = 0). Scalar
%                theta.corrinfo.corr = Correlation indicator. String. "uncorr" or "corr"
%
% Output:
%   meanV -- mean voltage (mV)
%   varV -- variance of voltage (mV^2)
%   outinfo -- extra info (struct)
%


%Load in parameters
Ke = theta.K.ee; %# Excitatory Synapses
Ki = theta.K.ei; %# Inhibitory Synapses
lam_e = theta.r.ee; % Excitatory Rate
lam_i = theta.r.ei; % Inhibitory Rate
W = theta.Wf; We = W(1,1); Wi = W(1,2); %Synaptic Weights, excitatory and inhibitory respectively
Ve = theta.ereversal; % Excitatory Reversal
Vi = theta.ireversal; % Inhibitory Reversal
tauleak = theta.tauleak; % Leak Time Constant
if isfield(theta,'I'); I = theta.I; else; I = 0; end %External Current
if isfield(theta.corrinfo,'corridx'); arho = theta.corrinfo.corridx; else; arho = 0; end
if ~isfield(theta,'ir'); irflag = false; else; irflag = theta.ir; end %Default no intrinisc resistance

%If the rate was inhomogenous, take the average
if length(lam_e) > 1; lam_e = mean(lam_e(:)); end
if length(lam_i) > 1; lam_i = mean(lam_i(:)); end

% Condition set for input type (Allows computation to be quicker for E and
% I only cases
if strcmp(theta.corrinfo.corr,'uncorr'); inputcase = "UnCorr"; alpha_e = 0; alpha_i = 0; beta_e = 0; beta_i = 0;
else
    alpha_e = theta.corrinfo.ae1; beta_e = theta.corrinfo.ae2;  %Excitatory Beta() Paramters
    alpha_i = theta.corrinfo.ai1; beta_i = theta.corrinfo.ai2;  %Inhibitory Beta() Paramters
    if Ki == 0 || beta_i == 0; inputcase = "Eonly";
    elseif Ke == 0 || beta_e == 0; inputcase = "Ionly";
    else; inputcase = "EI";
    end
end

%Adjusted Rates (b) from Raw: b = r*beta*(psi(beta+K)-psi(beta))
if alpha_e > 0
    ftemp = @(k) k.*exp((gammaln(Ke + 1)-gammaln(k + 1)-gammaln(Ke - k + 1)) +...
                           log(beta(k+alpha_e,(beta_e + Ke - k))./beta(alpha_e,beta_e)));
    eftemp = sum(ftemp(1:Ke),'omitnan'); 
    be = Ke*lam_e*1e-3/(eftemp);
else
    be = lam_e*beta_e*(psi(beta_e+Ke)-psi(beta_e)).*1e-3;
end

if alpha_i > 0
    ftemp = @(k) k.*exp((gammaln(Ki + 1)-gammaln(k + 1)-gammaln(Ki - k + 1)) +...
                           log(beta(k+alpha_i,(beta_i + Ki - k))./beta(alpha_i,beta_i)));
    eftemp = sum(ftemp(1:Ki),'omitnan'); 
    bi = Ki*lam_i*1e-3/(eftemp);
else
    bi = lam_i*beta_i*(psi(beta_i+Ki)-psi(beta_i)).*1e-3;
end

%Total Rate seperated based on across correlations: b = a*max(be,bi)+(1-a)*(be+bi)
bmax = max(be,bi);
bei = lam_e*beta_e*(psi(beta_e+Ke+Ki)-psi(beta_e))*1e-3;
b = arho*bei + (1-arho).*(sum([be,bi],'omitnan'));


% Adjust based on size of input (Using bionmial coef approximation begins
% to differ at around K = 1500. At this point, switch to Symbolic
% calculations. NOTE SYMBOLIC IS MUCH SLOWER
%if Ke >= 1500 || Ki >= 1500; inputsize = "Large"; warning('Large K - Using Symbolic Calculations. Will Be Slow');
%else; inputsize = "Small"; end

%Set conditions based on a -> 0 or not. When a = 0 replace denomentor by:
%den = psi(beta + K) - psi(beta) instead of beta(alpha,beta)
if alpha_e == 0; den_e = psi(beta_e+Ke)-psi(beta_e); else; den_e = beta(alpha_e,beta_e); end
if alpha_i == 0; den_i = psi(beta_i+Ki)-psi(beta_i); else; den_i = beta(alpha_i,beta_i); end
if alpha_i == 0 && alpha_e == 0; den_ei = psi(beta_e+Ke+Ki)-psi(beta_e); else; den_ei = beta(alpha_e,beta_e); end

%If we want intrinsic voltage dependent conductance. Adjust tau and W
if irflag
    if ~isfield(theta,'meanv'); error("Need Desired Mean V"); end
    mv = theta.meanv;
    gv = theta.irf;
    gv = 1+gv(mv);
else
    gv = 1;
end

%Calculate Mean and Var Seperated by Case and Size of input
switch inputcase
    case "UnCorr"
        ae1 = lam_e*1e-3*Ke*tauleak*(1-exp(-We/tauleak));
        ai1 = lam_i*1e-3*Ki*tauleak*(1-exp(-Wi/tauleak));
        meanV = (Ve*ae1 + Vi*ai1 + I)/(gv + ae1 + ai1);
        %meanV = (lam_e*1e-3*Ke*We*Ve + Ie)/(gv + lam_e*1e-3*Ke*We);
        %meanV = (lam_e*1e-3*Ke*We*Ve)./(gv.*(gv + lam_e*1e-3*Ke*We));
        
        alphae = lam_e*1e-3*Ke*tauleak*(1-exp(-We/tauleak)).^2/2;
        alphai = lam_i*1e-3*Ki*tauleak*(1-exp(-Wi/tauleak)).^2/2;
        
        varV = (alphae*(Ve-meanV).^2 + alphai*(Vi-meanV).^2)./(gv + ae1 + ai1);
        %varV = (alphae*((1+ai1)*Ve - ai1*Vi - Ie).^2 + alphai*((1+ae1)*Vi - ae1*Ve - Ie)^2)/...
        %        ((1 + ae1 + ai1).^2*(1 + gv + ae1 + ai1 - alphae - alphai));
        cei = 0;
    case "Eonly" %Excitatory Only Input Case
        epdf = @(k) exp(gammaln(Ke+1) + gammaln(beta_e+Ke-k) - gammaln(beta_e+Ke) - gammaln(Ke-k+1)-log(k))./den_e;
        
        %Marginal over E for [1-e^(We*k)/tau]
        fe1 = @(k) (1-exp(-(We.*k)./tauleak)).*epdf(k);
        efe1 = sum(fe1(1:Ke),'omitnan');

        %Marignal over E for [(1-e^(We*k)/tau)^2]
        fe2 = @(k) (1-exp(-((We.*k)./tauleak))).^2.*epdf(k);
        efe2 = sum(fe2(1:Ke),'omitnan');
        
        ae1 = b*tauleak*efe1; %a_{e1}
        meanV = (ae1*Ve+I)/(1 + ae1); %Mean V

        alphae = b*tauleak/2*efe2; %alpha_{e1}
        varV = (alphae*(Ve-meanV)^2)/(1 + ae1 - alphae); %Var V   
        
        ai1 = 0; alphai = 0; cei = 0; %No Inhibitory Components
    case "Ionly" %Inbhitory Only Input Case
        ipdf = @(k) exp(gammaln(Ki+1) + gammaln(beta_i+Ki-k) - gammaln(beta_i+Ki) - gammaln(Ki-k+1)-log(k))./den_i;
        
        %Marginal over I for [1-e^(Wi*k)/tau]
        fi1 = @(k) (1-exp(-(Wi.*k)./tauleak)).*ipdf(k);
        efi1 = sum(fi1(0:Ki),'omitnan');

        %Marignal over I for [(1-e^(Wi*k)/tau)^2]
        fi2 = @(k) (1-exp(-((Wi.*k)./tauleak))).^2.*ipdf(k);
        efi2 = sum(fi2(0:Ki),'omitnan');
        
        ai1 = b*tauleak*efi1; %a_{i1}
        meanV = (ai1*Vi)/(1+ai1); %Mean V

        alphai = b*tauleak/2*efi2; %alpha_{i1}
        varV = (alphai*(Vi-meanV)^2)/(1 + ai1 - alphai); %Var V
        
        ae1 = 0; alphae = 0; cei = 0; %No Excitatory Components
    case "EI" %Both Excitatory and Inhibitory Input
        f_coef1 = @(k1,k2,w1,w2) (k1.*w1./(k1.*w1+k2.*w2)).*(1-exp(-(w1.*k1 + w2.*k2)/tauleak)); %[k1/(k1+k2) * (1-exp(-(w1*k1 + w2*k2)/tau))]
        f_coef2 = @(k1,k2,w1,w2) (k1.*w1./(k1.*w1+k2.*w2)).*(1-exp(-(w1.*k1 + w2.*k2)/tauleak)).^2; %[k1/(k1+k2) * (1-exp(-(w1*k1 + w2*k2)/tau))^2]
        f_coefbeta = @(k1,k2,w1,w2) ((k1.*w1.*k2.*w2)./(k1.*w1 + k2.*w2).^2).*(1-exp(-(w1.*k1 + w2.*k2)/tauleak)).^2; %[k1*k2/(k1+k2)^2 * (1-exp(-(w1*k1 + w2*k2)/tau))]

        bnorme = be/(be+bi); bnormi = bi/(be+bi); %Coefficients for marginals
        
        coef1 = arho.*bmax/((1-arho).*(be+bi) + arho.*bmax); %Addjusted Coefficient Depending on Within Correlation (= 0 if a = 0)
        coef2 = (1-arho).*(be+bi)/((1-arho).*(be+bi) + arho.*bmax); %(=0 if a = 1)

        [WWE, WWI] = meshgrid(0:Ke,0:Ki);    
        fE_val = zeros(size(WWE)); fI_val = zeros(size(WWI)); 
        
        
        epdf = @(w) exp(gammaln(Ke+1) + gammaln(beta_e+Ke-w) - gammaln(beta_e+Ke) - gammaln(Ke-w+1)-log(w))./den_e;
        ipdf = @(w) exp(gammaln(Ki+1) + gammaln(beta_i+Ki-w) - gammaln(beta_i+Ki) - gammaln(Ki-w+1)-log(w))./den_i;
        eipdf = @(we,wi) exp(gammaln(Ke+1) + gammaln(Ki+1) + gammaln(alpha_e + we + wi) + gammaln(beta_e+Ke+Ki-we-wi)...
                                               -gammaln(we+1)-gammaln(Ke-we+1)-gammaln(wi+1)-gammaln(Ki-wi+1)-gammaln(alpha_e + beta_e + Ke + Ki))./den_ei;
        if arho == 0
            %Marginal over E for [1-e^(We*k)/tau]
            fE1_val = bnorme.*epdf(0:Ke); fE1_val(1,1) = 0; 
            efe1 = sum(f_coef1(0:Ke,0,We,Wi).*fE1_val,'all','omitnan');

            %Marginal over I for [1-e^(Wi*k)/tau]
            fI1_val = bnormi.*ipdf(0:Ki); fI1_val(1,1) = 0;
            efi1 = sum(f_coef1(0:Ki,0,Wi,Wi).*fI1_val,'all','omitnan');

            %Marignal over E for [(1-e^(We*k)/tau)^2]
            efe2 = sum(f_coef2(0:Ke,0,We,Wi).*fE1_val,'all','omitnan');

            %Marignal over I for [(1-e^(Wi*k)/tau)^2]
            efi2 = sum(f_coef2(0:Ki,0,Wi,We).*fI1_val,'all','omitnan');
        
            efcei = 0; %NO CROSS, CROSS TERM = 0

        else
            fE1_val = zeros(size(WWE)); fI1_val = zeros(size(WWI)); 

            %Marginal over E for [1-e^(We*k)/tau]
            fE1_val(1,:) = bnorme.*epdf(0:Ke); 
            %Marginal over I for [1-e^(Wi*k)/tau]
            fI1_val(:,1) = bnormi.*ipdf(0:Ki); 
           
            %Ratio between margianls (depends on who has the larger
            %rate, e or i)
            if be >= bi; bc1 = (be-bi)/be; bc2 = bi/be; fm_val = fE1_val/bnorme;
            elseif bi > be; bc1 = (bi-be)/bi; bc2 = be/bi; fm_val = fI1_val/bnormi;
            end
                    
            %fEI(ne_,ke,ni_,ki) = nchoosek(ne_,ke).*nchoosek(ni_,ki).*beta(alpha_e+ke+ki,beta_e+ne_+ni_-ke-ki)./den_ei;
            fEI_val = eipdf(WWE,WWI); fEI_val(1,1) = 0;
            fEI_val(isnan(fEI_val)) = 0;
            
            rhs = coef1.*(bc1.*fm_val + bc2.*fEI_val) + coef2.*(fE_val + fI_val);
            rhs(isinf(rhs)) = 0;
            efe1 = sum(f_coef1(WWE,WWI,We,Wi).*rhs,'all','omitnan');
            efi1 = sum(f_coef1(WWI,WWE,Wi,Wi).*rhs,'all','omitnan');

            efe2 = sum(f_coef2(WWE,WWI,We,Wi).*rhs,'all','omitnan');
            efi2 = sum(f_coef2(WWI,WWE,Wi,We).*rhs,'all','omitnan');
            efcei = sum(f_coefbeta(WWE,WWI,We,Wi).*rhs,'all','omitnan');
        end
        %}
        ae1 = b*tauleak*efe1;
        ai1 = b*tauleak*efi1;
        meanV = (Ve*ae1 + Vi*ai1+I)/(gv + ae1 + ai1);
        
        alphae = b*tauleak/2*efe2; alphai = b*tauleak/2*efi2;

        cei = b*tauleak/2*efcei;
        varV = (alphae*(Ve-meanV)^2 + alphai*(Vi-meanV)^2 - cei*(Ve-Vi).^2)/(gv + ae1 + ai1 - alphae - alphai);   
end

outinfo.ae1 = ae1;
outinfo.ai1 = ai1;
outinfo.alphae1 = alphae;
outinfo.alphai1 = alphai;
outinfo.cei = cei;
outinfo.tautest = tauleak;
outinfo.Wtest = We;
outinfo.Itest = I;
