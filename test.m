clear;

N = 10000;
E0 = [50 100:100:500 750 1000];
e = {};

for i = 1:length(E0)
    e{i} = simulateSEE(N,E0(i),'Si_DFT');
end

%% Histograms
% spec_se = zeros(1);
% spec_pe = zeros(1);
sey = zeros(size(E0));
bse = zeros(size(E0));
for i = 1:length(E0)
    bse(i) = 0;
    sey(i) = 0;
    for j = 1:length(e{i})
        if ~e{i}(j).Inside && ~e{i}(j).Dead
            if e{i}(j).isSecondary
                % if length(spec_se) == 1 && n_se == 0
                    % spec_se(1) = e(i).Energy;
                %     n_se = n_se + 1;
                % else
                    % spec_se(end+1) = e(i).Energy;
                    sey(i) = sey(i) + 1;
                % end
            else
                % if length(spec_se) == 1 && n_pe == 0
                %     spec_pe(1) = e(i).Energy;
                    % n_pe = n_pe + 1;
                % else
                    % spec_pe(end+1) = e(i).Energy;
                    bse(i) = bse(i) + 1;
                % end
            end
        end
    end
end

%% Yields
% {
figure
hold on
load('C:\Users\onr5\OneDrive - NIST\dev\m-scripts\jmonsel_si.mat')
plot(jmonsel(:,1),jmonsel(:,2),DisplayName='JMONSEL')
plot(E0,(sey+bse)/N,DisplayName='TEY')
plot(E0,bse/N,DisplayName='BSE')
plot(E0,sey/N,DisplayName='SEY')
legend
%}

%% Plots
%{
figure
hold on
histogram(spec_pe,500)
histogram(spec_se,500)
xlabel('Electron energy (eV)')
%}

