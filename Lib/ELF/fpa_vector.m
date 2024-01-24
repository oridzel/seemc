function [elf,elf_pl,elf_se] = fpa_vector(q,omega,optical_omega,optical_elf)

q(isnan(q)) = 0;
omega_pl = eps:0.01:1e4;
omega_pl = omega_pl/h2ev;
elf_pl = zeros(length(omega),size(q,2));
elf_se = zeros(length(omega),size(q,2));
omega_0 = zeros(size(omega));

for i = 1:size(q,2)
    epsilon = lindhard(q(:,i),omega,omega_pl);
    q_m = q_minus(omega',omega_pl);
    q_p = q_plus(omega',omega_pl);
    g = g_coef(omega_pl,optical_omega,optical_elf);
    im = imag(epsilon)./(imag(epsilon).^2 + real(epsilon).^2);
    se = g .* im .* heaviside(q_p - q(:,i)) .* heaviside(q(:,i) - q_m);
    se(isnan(se)) = 0;

    for j = 1:length(omega)
        if size(q,1) > 1
            guess_interval = [0,min(omega(j),abs(q(j,i)/2-omega(j)/q(j,i))) + 1e-5];
            fun = @(x) real_lindhard(q(j,i),omega(j),x);
        else
            guess_interval = [0,min(omega(j),abs(q(i)/2-omega(j)/q(i))) + 1e-5];
            fun = @(x) real_lindhard(q(i),omega(j),x);
        end
        try
            omega_0(j) = fzero(fun,guess_interval);
        catch
            omega_0(j) = 0;
        end        
    end

    elf_pl(:,i) = g_coef(omega_0,optical_omega,optical_elf)' .* pi./abs(d_epsilon_real(q(:,i),omega',omega_0')) .* heaviside(q_minus(omega',omega_0') - q(:,i)); 
    elf_pl(isnan(elf_pl)) = 0;
    
    elf_se(:,i) = trapz(omega_pl,se,2);
    elf = elf_pl + elf_se;
end

    function val = real_lindhard(q,omega,omega_pl)
        epsilon_lindhard = lindhard(q,omega,omega_pl);
        val = real(epsilon_lindhard);
    end
    
    function val = d_epsilon_real(q,omega,omega_pl)
        kf = k_f(omega_pl);
        x = 2*omega ./ k_f(omega_pl).^2;
        z = q ./ (2*kf);
        u = x ./ (4*z);
        y_plus = z + u;
        y_minus = z - u;
    
        val = ( log(abs( (y_minus + 1)./(y_minus - 1)) ) + log(abs( (y_plus + 1)./(y_plus - 1)) ) ) ./ (4*kf.^(5/2)*sqrt(3*pi).*z.^3);
        ind_1 = x > 100*z;
        ind_2 = z > 100*x;
    
        if any(unique(ind_1))
            a = z ./ x;
            d = 16 ./ (kf.^(5/2)*sqrt(3*pi).*x.^2) .* (-1 - 16.*a.^2 - 16.*a.^4.*(16 + x.^2) - 512/3*a.^6.*(24 + 5*x.^2));
            val(ind_1) = d(ind_1);
        end
        if any(unique(ind_2))
            b = x ./ (z.*(z.^2 - 1));
            d = ( log( ((z+1)./(z-1)).^2 ) + 4.*z.*b.^2.*( 1 + (1 + z.^2).*b.^2 + 1/3*(3 + z.^2).*(1 + 3*z.^2).*b.^4 ) ) ./ (4*kf.^(5/2)*sqrt(3*pi).*z.^3);
            val(ind_2) = d(ind_2);
        end
    end
    
    function epsilon = lindhard(q,omega,omega_pl)
        kf = k_f(omega_pl);
        x = 2*omega' ./ kf.^2;
        z = q ./ (2*kf);
        u = x ./ (4*z);
    
        if all(all(isnan(u)))
            epsilon_real = ones(size(omega_pl));
            epsilon_imag = zeros(size(omega_pl));
        elseif all(all(z == 0)) && all(all(x ~= 0))
            epsilon_real = 1 - 16 ./ (3*kf*pi.*x.^2);
            epsilon_imag = zeros(size(epsilon_real));
        else
            epsilon_real = 1 + 1./(pi.*kf.*z.^2) .* ( 1/2 + 1./(8*z).*(f(z - u) + f(z + u)) );
            ind = z == 0;
            e_1 = 1 - 16 ./ (3*kf*pi.*x.^2);
            epsilon_real(ind) = e_1(ind);
            coef = 1 ./ (8*kf.*z.^3);
            epsilon_imag = zeros(size(epsilon_real));
            ind_1 = x > 0 & x < 4*z.*(1 - z);
            ind_2 = x > abs(4*z.*(1 - z)) & x < 4*z.*(1 + z);
            if any(unique(ind_1))
                e_2 = coef.*x;
                epsilon_imag(ind_1) = e_2(ind_1);
            end
            if any(unique(ind_2))
                e_2 = coef.*(1 - (z - u).^2);
                epsilon_imag(ind_2) = e_2(ind_2);
            end
        
            ind_1 = u < 0.01;
            u_over_z = u./(z+1);
            ind_2 = u_over_z > 100;
        
            if any(unique(ind_1))
                e_1 = 1 + 1./(pi.*kf.*z.^2) .* ( 1/2 + 1./(4*z).*( (1 - z.^2 - u.^2).*log(abs((z+1)./(z-1))) + (z.^2 - u.^2 - 1)*2.*u.^2.*z./((z.^2 - 1).^2) ) );
                e_2 = u./(2*kf.*z.^2);
                epsilon_real(ind_1) = e_1(ind_1);
                epsilon_imag(ind_1) = e_2(ind_1);
            end
        
            if any(unique(ind_2))
                e_1 = 1 - 16 ./ (3*kf*pi.*x.^2) - 256*z.^2 ./ (5*kf*pi.*x.^4) - 256*z.^4 ./ (3*kf*pi.*x.^4);
                epsilon_real(ind_2) = e_1(ind_2);
                epsilon_imag(ind_2) = 0;
            end
        end
        
        epsilon = complex(epsilon_real,epsilon_imag);
    end
    
    function val = f(t)
        val = (1 - t.^2).*log(abs((1+t)./(-1+t)));
        val(abs(t) == 1) = 0;
    end
    
    function val = q_plus(omega,omega_pl)
        val = k_f(omega_pl) + sqrt(k_f(omega_pl).^2 + 2*omega);
    end
    
    function val = q_minus(omega,omega_pl)
        val = (-1)*k_f(omega_pl) + sqrt(k_f(omega_pl).^2 + 2*omega);
    end
    
    function val = k_f(omega)
        val = (3*pi/4)^(1/3) * omega.^(2/3);
    end
    
    function g = g_coef(omega_pl,optical_omega,optical_elf)  
        g = 2./(pi*omega_pl) .* interp1(optical_omega,optical_elf,omega_pl*h2ev);
        g(omega_pl < optical_omega(1)) = eps;
    end

end