clear;

ef = 12;
de = 30;

for i = 1:10000
    e(i) = fegdos(de,ef);
end

histogram(e)