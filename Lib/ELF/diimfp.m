function [iimfp,x_in] = diimfp(E0,eloss,optical_eloss,optical_elf,is_metal,egap,varargin)
    if nargin < 5
        is_metal = true;
        egap = 0;
    end

    if is_metal
        E0 = E0/h2ev;
    else
        E0 = (E0 - egap)/h2ev;
    end
    omega = eloss/h2ev;
    n_q = 100;   
    C = 137.036; % a.u.
    q1 = log( sqrt(E0*(2 + E0 / C^2 )) - sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q2 = log( sqrt(E0*(2 + E0 / C^2 )) + sqrt((E0 - omega).*(2.0 + (E0 - omega) / C^2 )) );
    q = zeros(length(omega),n_q);   
    for i = 1:n_q
        q(:,i) = q1 + (i-1)*(q2-q1)/n_q;
    end
    res = fpa_vector(exp(q),omega,optical_eloss,optical_elf);
    
    iimfp = zeros(size(eloss));
    iimfp(1) = eps;
    for i = 2:length(eloss)
        iimfp(i) = trapz(q(i,:),res(i,:));
    end
    
    rel_cor_factor = (1 + E0/C^2)^2 / (1 + E0/(2*C^2));
    iimfp = rel_cor_factor / (pi*E0) * iimfp;
    x_in = iimfp / (h2ev * a0); % now in 1/(eV A)
end