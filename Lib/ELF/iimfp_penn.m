function iimfp = iimfp_penn(E0, eloss, ELF)
%%
%{
   Calculates the IMFP values from optical data.
   The Penn algoritm (SSPA) is used.
   \param [in] E0 - the energy with respect to the Fermi level for which the diimfp is to be calculated
%}
%%

w = eloss/h2ev;

x_in = zeros(length(w),1);
    
    for i = 1:length(w)
        if w(i)<0.05
            x_in(i) = 0;
        else
            hhw=w(i);
            wpmax=2*(hhw-E0/h2ev+sqrt(E0/h2ev*(E0/h2ev-hhw)));
            i1=0;
            k=2;
            quit=false;
            while ~quit
                if w(k)>wpmax
                    quit=true;
                    wup=wpmax;
                else
                    wup=w(k);
                end
                wlo=w(k-1);
                if wlo>wpmax
                    wlo=wpmax;
                end
                temp = ELF(k)*(ddmfpp_integrand(hhw,wup)-ddmfpp_integrand(hhw,wlo));
                if isnan(temp)
                    temp = 0;
                end
                i1 = i1 + temp;
                k = k+1;
                x_in(i) = i1/2.0/pi/hhw/a0/E0;
            end
        end
    end
    iimfp = 1/trapz(eloss,x_in);
end

function x = ddmfpp_integrand(de,hwp)
    if hwp>de
        error('ddmfpp_integrand:Somethings wrong');
    else
        x = -(hwp+de*log(de-hwp));
    end
end



