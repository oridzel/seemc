function W = Make_W
%% W
W.Mat = 'W';
W.M = 183.84;
W.Z = 74;
W.Density = 19.25; %g/cm^3
W.Density = W.Density*10^-21/W.M*6.022*10^23;
W.NvTPP = 6;
W.NvSGS = 6;
W.Eg = 0;
W.Ep = 22.86; %free-electron plasmon energy
W.Ef = 10.1; %Fermi energy

E0 = [100:20:500 600:100:2500 2750:250:5000 5500:500:40000]';
W.Elastic.x = zeros(numel(E0),1);
W.Elastic.l_el = zeros(numel(E0),1);
W.Elastic.l_tr = zeros(numel(E0),1);
W.Elastic.x = E0;
W.DECS.E0 = E0;

tic;
[data] = ElsepaRunner.RunElsepa(W.Z,E0);
toc
W.DECS.x = data(1).x;
for i = 1:numel(E0)
    W.Elastic.l_el(i) = 1/data(i).sigma_el/W.Density;
    W.Elastic.l_tr(i) = 1/data(i).sigma_tr1/W.Density;
    W.DECS.y(:,i) = data(i).y/trapz(data(i).x,data(i).y);
end

W.DIIMFP.E0 = [100;200;500;1000;2000;5000;10000;];
for i=1:numel(W.DIIMFP.E0)
    WernerData = load([cd '/W_in/' W.Mat num2str(W.DIIMFP.E0(i)) '.diimfp']);
    W.DIIMFP.y(:,i) = WernerData(:,3);
    if i==1
         W.DIIMFP.x = WernerData(:,1);
    end
end

W.XPS.E = [100; 200; 500; 1000; 1500; 2000; 3000; 4000; 5000];
W.XPS.Shells = {'3D3/2','3D5/2', '4S1/2', '4P1/2', '4P3/2', '4D3/2', '4D5/2', '4F5/2', '4F7/2', '5S1/2', '5P1/2', '5P3/2', '5D3/2', '6S1/2'};
W.XPS.EB = [1874; 1812; 599;495;428;261;248;39;36;80;51;41;9;8];
W.XPS.PCS = [0.3592e3 0.3239e3 0.230e3 0.1389e3 0.8903e2 0.5991e2 0.3029e2 0.1706e2 0.1038e2;...
    0.5577e3 0.4995e3 0.3504e3 0.2066e3 0.1304e3 0.8668e2 0.4297e2 0.2384e2 0.1432e2;
    0.6184e2 0.5389e2 0.3706e2 0.2230e2 0.1487e2 0.1061e2 0.6176e1 0.4025e1 0.2823e1;
    0.6021e2 0.5639e2 0.4393e2 0.2856e2 0.1952e2 0.1400e2 0.8028e1 0.5090e1 0.3458e1;
    0.2371e3 0.2034e3 0.1328e3 0.7395e2 0.4625e2 0.3124e2 0.1653e2 0.9925e1 0.6470e1;
    0.2359e3 0.2886e3 0.2333e3 0.1182e3 0.6398e2 0.3770e2 0.1579e2 0.7852e1 0.4374e1;
    0.3894e3 0.4680e3 0.3600e3 0.1752e3 0.9277e2 0.5381e2 0.2206e2 0.1080e2 0.5941e1;
    0.3590e4 0.2986e4 0.8356e3 0.1584e3 0.4828e2 0.1916e2 0.4707e1 0.1626e1 0.6889e0;
    0.4965e4 0.4011e4 0.1092e4 0.2037e3 0.6156e2 0.2429e2 0.5911e1 0.2029e1 0.8546e0;
    0.7806e2 0.5059e2 0.2040e2 0.8372e1 0.4660e1 0.2991e1 0.1542e1 0.9387e0 0.6301e0;
    0.7009e2 0.3397e2 0.1604e2 0.7802e1 0.4652e1 0.3077e1 0.1606e1 0.9648e0 0.6329e0;
    0.1928e3 0.1127e3 0.4580e2 0.1884e2 0.1028e2 0.6412e1 0.3099e1 0.1769e1 0.1116e1;
    0.1356e3 0.2470e2 0.2513e2 0.1169e2 0.5911e1 0.3329e1 0.1319e1 0.6344e0 0.3458e0;
    0.1767e2 0.8361e1 0.2557e1 0.9217e0 0.4894e0 0.3068e0 0.1543e0 0.9277e-1 0.6179e-1]'*10^-7; %kbarn
W.XPS.betta = [0.371 0.462 0.771 1.031 1.148 1.204 1.237 1.223 1.190;
    0.424 0.551 0.859 1.082 1.170 1.204 1.205 1.171 1.126;
    1.812 1.834 1.867 1.889 1.899 1.905 1.915 1.923 1.930;
    0.101 0.627 1.239 1.527 1.621 1.662 1.683 1.675 1.655;
    0.143 0.618 1.187 1.493 1.613 1.676 1.731 1.746 1.744;
    -0.801 -0.206 0.701 1.143 1.293 1.351 1.369 1.337 1.288;
    -0.697 0.001 0.860 1.223 1.325 1.351 1.328 1.272 1.210;
    0.390 0.424 0.777 1.003 1.048 1.034 0.953 0.863 0.780;
    0.381 0.435 0.794 1.009 1.045 1.026 0.941 0.853 0.772;
    1.675 1.781 1.860 1.890 1.901 1.908 1.916 1.924 1.930;
    0.057 0.160 1.159 1.554 1.658 1.696 1.711 1.698 1.674;
    -0.258 0.325 1.140 1.504 1.631 1.695 1.748 1.761 1.756;
    1.977 -0.680 0.615 1.167 1.325 1.383 1.395 1.358 1.306;
    1.623 1.764 1.858 1.889 1.901 1.908 1.916 1.923 1.930;]';
W.XPS.gamma = [-0.676e-1 -0.124e0 -0.215e0 -0.152e0 -0.812e-2 0.148e0 0.443e0 0.694e0 0.910e0;
    -0.683e-1 -0.129e0 -0.216e0 -0.129e0 0.284e-1 0.190e0 0.484e0 0.727e0 0.930e0;
    0.622e0 0.632e0 0.597e0 0.489e0 0.378e0 0.275e0 0.108e0 -0.100e-1 -0.872e-1;
    -0.163e0 0.107e0 0.420e0 0.317e0 0.156e0 0.511e-1 -0.150e-1 0.297e-1 0.126e0;
    -0.119e0 0.896e-1 0.312e0 0.213e0 0.776e-1 -0.605e-2 -0.395e-1 0.312e-1 0.150e0;
    -0.151e0 -0.155e-1 0.133e0 -0.282e-1 -0.326e-1 0.472e-1 0.284e0 0.526e0 0.745e0;
    -0.147e0 0.366e-1 0.127e0 -0.365e-1 -0.216e-1 0.720e-1 0.319e0 0.558e0 0.769e0;
    0.153e-2 0.736e-3 0.245e-1 0.143e0 0.302e0 0.455e0 0.711e0 0.900e0 0.104e1;
    0.185e-2 0.139e-2 0.257e-1 0.148e0 0.311e0 0.467e0 0.723e0 0.910e0 0.105e1;
    0.170e0 0.303e0 0.421e0 0.408e0 0.338e0 0.256e0 0.106e0 -0.706e-2 -0.833e-1;
    0.158e0 -0.507e-1 0.315e0 0.314e0 0.170e0 0.622e-1 -0.173e-1 0.175e-1 0.108e0;
    -0.566e-2 -0.279e-1 0.246e0 0.218e0 0.931e-1 0.540e-2 -0.418e-1 0.189e-1 0.132e0;
    0.268e0 -0.112e0 0.177e0 0.622e-3 -0.339e-1 0.296e-1 0.255e0 0.497e0 0.717e0;
    0.107e0 0.247e0 0.401e0 0.405e0 0.337e0 0.256e0 0.108e0 -0.540e-2 -0.823e-1]';
W.XPS.delta = [-0.594e-2 -0.190e-1 -0.108e-1 0.258e-1 0.496e-1 0.658e-1 0.890e-1 0.107e0 0.124e0;
    -0.757e-2 -0.200e-1 -0.115e-1 0.227e-1 0.455e-1 0.618e-1 0.874e-1 0.109e0 0.130e0;
    0.210e-2 0.146e-2 0.415e-3 -0.390e-3 -0.851e-3 -0.118e-2 -0.164e-2 -0.195e-2 -0.218e-2;
    0.112e0 0.913e-1 0.356e-1 0.706e-2 -0.451e-3 -0.306e-2 -0.456e-2 -0.427e-2 -0.285e-2;
    0.957e-1 0.791e-1 0.318e-1 0.636e-2 0.140e-2 0.154e-2 0.478e-2 0.807e-2 0.110e-1;
    0.116e0 0.130e0 0.248e-1 0.475e-2 0.156e-1 0.274e-1 0.477e-1 0.655e-1 0.829e-1;
    0.113e0 0.114e0 0.208e-1 0.519e-2 0.148e-1 0.258e-1 0.468e-1 0.672e-1 0.881e-1;
    0.807e-2 0.154e-1 0.344e-1 0.642e-1 0.909e-1 0.115e0 0.160e0 0.203e0 0.245e0;
    0.698e-2 0.148e-1 0.345e-1 0.643e-1 0.908e-1 0.115e0 0.161e0 0.205e0 0.248e0;
    0.195e-2 0.129e-2 0.454e-3 -0.232e-3 -0.672e-3 -0.101e-2 -0.148e-2 -0.182e-2 -0.206e-2;
    -0.297e-2 0.379e-1 0.184e-1 0.821e-3 -0.304e-2 -0.423e-2 -0.512e-2 -0.508e-2 -0.408e-2;
    0.519e-2 0.365e-1 0.199e-1 0.222e-2 -0.951e-3 -0.182e-3 0.337e-2 0.661e-2 0.935e-2;
    0.397e-3 0.578e-1 0.179e-1 -0.325e-3 0.102e-1 0.217e-1 0.417e-1 0.594e-1 0.766e-1;
    0.169e-2 0.123e-2 0.486e-3 -0.200e-3 -0.647e-3 -0.980e-3 -0.146e-2 -0.180e-2 -0.204e-2;]';

