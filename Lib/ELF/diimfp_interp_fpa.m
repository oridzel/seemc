function [iimfp,x_in] = diimfp_interp_fpa(E0,eloss,fpa_elf,fpa_q,fpa_omega,osc)
    E0 = (E0 - osc.egap)/h2ev;
    omega = eloss/h2ev;
    n_q = 100;   
    C = 137.036; % a.u.
    q1 = log( sqrt(E0*(2 + E0 / C^2 )) - sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q1(E0 - omega == E0) = 0;
    q2 = log( sqrt(E0*(2 + E0 / C^2 )) + sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q = zeros(length(omega),n_q);   
    for i = 1:n_q
        q(:,i) = q1 + (i-1)*(q2-q1)/n_q;
    end

    [x,y] = meshgrid(fpa_q,fpa_omega);
    www = repmat(eloss,n_q,1);
    qqq = transpose(exp(q)/a0);
    qqq(isnan(qqq)) = 0;
    res_interp = transpose(interp2(x,y,fpa_elf,qqq,www));

    ind = isnan(res_interp);
    if any(unique(ind))      
        osc.qtran = exp(q)/a0;
        osc.eloss = eloss;
        res_mermin = eps_sum_allwq(osc,'bulk');
        res_interp(ind) = res_mermin(ind);
    end
    res = res_interp;
    
    iimfp = zeros(size(eloss));
    iimfp(1) = eps;
    for i = 2:length(eloss)
        iimfp(i) = trapz(q(i,:),res(i,:));
    end
    
    rel_cor_factor = (1 + E0/C^2)^2 / (1 + E0/(2*C^2));
    iimfp = rel_cor_factor / (pi*E0) * iimfp;
    x_in = iimfp / (h2ev * a0); % now in 1/(eV A)
end