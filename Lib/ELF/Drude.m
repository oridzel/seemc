function [eps_re, eps_im]=Drude(q,w_global,omega0,gamma,alpha,FermiEnergy)

    w_global = w_global(:);

    if alpha < 100
        % now there are different aplha's in use in the literature here we use dispersion relative to free particle
        %here not really important, dispersion is different than expected in Drude
        w_at_q = omega0 + 0.5*alpha * q.^2;
    else
        % now we use full dispersion if alpha > 100 with the Fermi energy based on the total electron density, so no alpha fitting possible
        w_at_q_square = omega0 * omega0 + (2.0 / 3.0)*FermiEnergy*q.^2 + q.^4/ 4.0;
        w_at_q = sqrt(w_at_q_square);
    end
    
    mm = w_global.^2 - w_at_q.^2;
    divisor = mm.^2 + w_global.^2*gamma^2;
  
    eps_re = mm ./ divisor;
    eps_im = w_global*gamma ./ divisor;
    
end