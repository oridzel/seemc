function base_uvw = updateDirection(uvw, deflection, algorithm)
    if algorithm == 1
        sinpsi = sin(deflection(1));
        cospsi = cos(deflection(1));
        sinfi = sin(deflection(2));
        cosfi = cos(deflection(2));
        costh = uvw(3);
        sinth = sqrt(uvw(1)^2 + uvw(2)^2);
        if sinth > 1e-10
            cosphi = uvw(1)/sinth;
            sinphi = uvw(2)/sinth;
        else
            cosphi = 1;
            sinphi = 0;
        end

        h0 = sinpsi*cosfi;
        h1 = sinth*cospsi + h0*costh;
        h2 = sinpsi*sinfi;
        base_uvw(1) = h1*cosphi - h2*sinphi;
        base_uvw(2) = h1*sinphi + h2*cosphi;
        base_uvw(3) = costh*cospsi - h0*sinth;
    elseif algorithm == 2
        theta0 = acos(uvw(end));
        phi0 = atan2(uvw(2),uvw(1));
        theta = acos(cos(theta0)*cos(deflection(1)) - sin(theta0)*sin(deflection(1))*cos(deflection(2)));
        phi = asin( sin(deflection(1))*sin(deflection(2))/sin(theta) ) + phi0;

        base_uvw(1) = sin(theta)*cos(phi);
        base_uvw(2) = sin(theta)*sin(phi);
        base_uvw(3) = cos(theta);
    end
end