function dsep = dsep_Tung(osc,E0,theta,decdigs,varargin)
%%
%{
   Calculates the normalised DSEP
   for a given energy and angle
   according to the Tung algorithm Eq.(9)
   C.J. Tung et al. / Phys. Rev. B 49 (1994) 16684.
%}

if nargin<4, decdigs=10; end
if nargin<3
    warning ('Error in input format')
else
    
    ind = bsxfun(@gt,osc.eloss,0);
    eloss = osc.eloss; 
    qmin = sqrt(2*E0/h2ev)-sqrt(2*(E0/h2ev-osc.eloss(ind)/h2ev));
    qmax = sqrt(2*E0/h2ev)+sqrt(2*(E0/h2ev-osc.eloss(ind)/h2ev));
    
    q = zeros(length(osc.eloss(ind)),2^(decdigs-1)+1);
    x_in_plus = zeros(size(osc.eloss));
    x_in_minus = zeros(size(osc.eloss));
    
    for i = 1:2^(decdigs-1)+1
        q(:,i) = qmin + (i-1)*(qmax-qmin)/2.^(decdigs-1);
    end
    
    if strcmp( osc.model,'Mermin')
        q(q==0) = 0.01;
    end
    
    osc.qtran = q/a0;
    osc.eloss = eloss(ind);
    
%     ELF = eps_sum_surf(osc);
    ELF = eps_sum_allwq(osc,'surface');
    
    q_s_plus = sqrt(q.^2 - ( (q.^2./2 + osc.eloss/h2ev) ./ sqrt(2*E0/h2ev) ).^2).*cosd(theta) + ( (q.^2./2 + osc.eloss/h2ev) ./sqrt(2*E0/h2ev)).*sind(theta); 
    q_s_minus = sqrt(q.^2 - ( (q.^2./2 + osc.eloss/h2ev) ./ sqrt(2*E0/h2ev) ).^2).*cosd(theta) - ( (q.^2./2 + osc.eloss/h2ev) ./sqrt(2*E0/h2ev)).*sind(theta);
    
    res_plus = bsxfun(@times,ELF,abs(q_s_plus))./(q.^3);
    res_minus = bsxfun(@times,ELF,abs(q_s_minus))./(q.^3);
    res_plus(isnan(res_plus))=0;
    res_minus(isnan(res_minus))=0;
    
    x_in_plus(1) = 0.0;
    x_in_minus(1) = 0.0;
    
    for i=2:length(osc.eloss)
        x_in_plus(i) = 1/pi/(E0/h2ev) * trapz(q(i,:),res_plus(i,:)) *(1/h2ev/a0);
        x_in_minus(i) = 1/pi/(E0/h2ev) * trapz(q(i,:),res_minus(i,:)) *(1/h2ev/a0);
    end

    dsep = x_in_plus + x_in_minus; %./trapz(osc.eloss,x_in);
end

end










