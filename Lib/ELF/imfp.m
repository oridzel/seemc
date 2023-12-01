function lambda_in = imfp(osc,energy,material,varargin)
    if nargin < 3
        material = 'metal';
    end

    if strcmp(material,'metal')
        osc.eloss = eps:0.1:energy - osc.Ef;
    elseif strcmp(material,'insulator')
        osc.eloss = osc.egap:0.1:energy - (osc.egap + osc.vb);
    end
    
    [iimfp,diimfp] = ndiimfp(osc,energy,false);
    lambda_in = 1/trapz(osc.eloss/h2ev,iimfp)*a0;
end