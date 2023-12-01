clear

dos = importdata('dos.in');
evb = dos(76,1);

d_ev = evb/200;
de = dos(:,1);

for i = 1:numel(de)
    for j = 1:200
        res(j,1) = j*d_ev;
        res(j,2) = interp1(dos(:,1),dos(:,2),de(i)+j*d_ev)*interp1(dos(:,1),dos(:,2),j*d_ev);
    end
    r(i,:,:) = [res(:,1),cumtrapz(res(:,1),res(:,2))];
end


e_int = interp2(res(:,1),de,r(:,:,2),res(:,1),20);