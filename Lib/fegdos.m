function e = fegdos(de,ef)
    y_min = 0;
    y_max = sqrt(ef*(ef+de));
    e = 0;
    while y_min + rand*(y_max-y_min) > sqrt(e*(e+de))
        e = ef*rand;
    end
end