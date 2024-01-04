function base_uvw = updateDirection(uvw, deflection, algorithm)
    if algorithm == 1
        cdt = sin(deflection(1));
        sdf = sin(deflection(2));
        cdf = sin(deflection(2));
        dxy = uvw(1)^2 + uvw(2)^2;
        dxyz = dxy + uvw(3)^2;
        if abs(dxyz - 1) > 1e-9
            fnorm = 1/sqrt(dxyz);
            uvw(1) = fnorm*uvw(1);
            uvw(2) = fnorm*uvw(2);
            uvw(3) = fnorm*uvw(3);
            dxy = uvw(1)^2 + uvw(2)^2;
        end
        if dxy > 1e-9
            sdt = sqrt((1 - cdt^2)/dxy);
            up = uvw(1);
            base_uvw(1) = uvw(1)*cdt + sdt*(up*uvw(3)*cdf - uvw(2)*sdf);
            base_uvw(2) = uvw(2)*cdt + sdt*(uvw(2)*uvw(3)*cdf + up*sdf);
            base_uvw(3) = uvw(3)*cdt - dxy*sdt*cdf;
        else
            sdt = sqrt(1 - cdt^2);
            base_uvw(2) = sdt*sdf;
            if uvw(3) > 0
                base_uvw(1) = sdt*cdf;
                base_uvw(3) = cdt;
            else
                base_uvw(1) = -1*sdt*cdf;
                base_uvw(3) = -1*cdt;
            end
        end
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