function [iimfp,diimfp] = ndiimfp(osc,E0,custom_elf,custom_q,custom_omega,varargin)

    if nargin < 3
        custom_elf = false;
    end

    E0 = (E0 - osc.egap)/h2ev;
    omega = osc.eloss/h2ev;
    n_q = 100;
    
    C = 137.036; % a.u.
    q1 = log( sqrt(E0*(2 + E0 / C^2 )) - sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q2 = log( sqrt(E0*(2 + E0 / C^2 )) + sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q = zeros(length(omega),n_q);
    
    for i = 1:n_q
        q(:,i) = q1 + (i-1)*(q2-q1)/n_q;
    end
    
    if strcmp( osc.model,'Mermin')
        q(q==0) = 0.01;
    end
    
    iimfp = zeros(size(osc.eloss));

    if custom_elf
        [x,y] = meshgrid(custom_q,custom_omega);
        www = repmat(osc.eloss,n_q,1);
        qqq = transpose(exp(q)/a0);
        res_custom = transpose(interp2(x,y,custom_elf,qqq,www));
        if E0*h2ev >= 110
            ind = isnan(res_custom);
            osc.qtran = exp(q)/a0;
            res_mermin = eps_sum_allwq(osc,'bulk');
            res_custom(ind) = res_mermin(ind);
        end
        res = res_custom;
    else
        osc.qtran = exp(q)/a0; 
        res = eps_sum_allwq(osc,'bulk');
    end

    rel_cor_factor = (1 + E0/C^2)^2 / (1 + E0/(2*C^2));
    
    iimfp(1) = eps;
    for i = 2:length(osc.eloss)
        iimfp(i) = trapz(q(i,:),res(i,:));
    end
    
    iimfp = rel_cor_factor*iimfp / (pi*E0);
    diimfp = iimfp./trapz(omega,iimfp) / (h2ev * a0); %normalized
end











