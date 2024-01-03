classdef ElsepaRunner
    
    methods (Static)
        function [Res] = RunElsepa(Composition, E0, custom_input_file, isMolecule, varargin)
            if nargin < 3
                custom_input_file = '';
                isMolecule = false;
            end
            if isempty(E0); return; end
            current_full_path = dbstack('-completenames');
            current_file_name = dbstack;
            ind = strfind(current_full_path(1).file,current_file_name(1).file);
            dir_Elsepa = [current_full_path(1).file(1:ind-2) filesep 'elsepa-2020'];
            if ~isempty(custom_input_file)
                filename = custom_input_file;
            else
                filename = 'lub.in';
            end
            f_name = fullfile(dir_Elsepa,filename);
            sumweights = sum(Composition.index);
            sigma_el = zeros(length(Composition.Z),numel(E0));
            sigma_tr1 = zeros(length(Composition.Z),numel(E0));

            for z = 1:length(Composition.Z)
            
                if isempty(custom_input_file)
                    fid = fopen(f_name,'w');
                    fprintf(fid, 'IZ      %1.0f         atomic number                               [none]\n',Composition.Z(z));
                    % fprintf(fid, 'MNUCL   3          rho_n (1=P, 2=U, 3=F, 4=Uu)                  [  3]\n');
                    % fprintf(fid, 'MELEC   4          rho_e (1=TFM, 2=TFD, 3=DHFS, 4=DF, 5=file)   [  4]\n');
                    fprintf(fid, 'MUFFIN  0          0=free atom, 1=muffin-tin model              [  0]\n');
                    % fprintf(fid, 'RMUF   -1        muffin-tin radius (cm)                  [measured]\n');
                    fprintf(fid, 'IELEC  -1          -1=electron, +1=positron                     [ -1]\n');
                    % fprintf(fid, 'MEXCH   1          V_ex (0=none, 1=FM, 2=TF, 3=RT)              [  1]\n');
    
                    for i=5:numel(E0)
                        fprintf(fid, 'EV      %1.0f       kinetic energy (eV)                         [none]\n',E0(i));
                    end
                    
                    fclose(fid);
                end

                cd(dir_Elsepa);
                if ispc
                    if isMolecule
                        system(['elsepm < ' filename]);
                    else
                        system(['elsepa < ' filename]);
                    end
                elseif ismac || isunix
                    if isMolecule
                        system(['./elsepm < ' filename]);
                    else
                        system(['./elsepa < ' filename]);
                    end
                end
                    
                DELIMITER = ' ';
            
                for i=1:numel(E0)
    
                    if E0(i) < 5
                        energy = 5;
                    else
                        energy = E0(i);
                    end
                    
                    f_name_el=fullfile(dir_Elsepa,...
                        ['dcs_' strrep(strrep(num2str(energy, '%1.3e'),'.','p'),'+0','0') '.dat']);
                    
                    data = struct;
    
                    if isMolecule
                        HEADERLINES = 34;
                    else
                        % For muffin-tin ON
                        % if E0(i) < 100
                        %     HEADERLINES = 47;
                        % else
                        %     HEADERLINES = 44;
                        % end
    
                        if E0(i) < 100
                            HEADERLINES = 39;
                        else
                            HEADERLINES = 36;
                        end
                    end
    
                    El = importdata(f_name_el, DELIMITER, HEADERLINES);
                    if isMolecule
                        a = El.textdata(~cellfun('isempty',regexp(El.textdata,'total cs')));
                        sigma_el(z,i) = str2double(extractBetween(a{1},": "," "))*1e16;
                        a = El.textdata(~cellfun('isempty',regexp(El.textdata,'1st transport cs')));
                        sigma_tr1(z,i) = str2double(extractBetween(a{1},": "," "))*1e16;
                    else
                        a = El.textdata(~cellfun('isempty',regexp(El.textdata,'Total elastic cross section')));
                        sigma_el(z,i) = str2double(extractBetween(a{1},"cm**2 = "," a0**2"))*a0^2;
                        a = El.textdata(~cellfun('isempty',regexp(El.textdata,'1st transport cross section')));
                        sigma_tr1(z,i) = str2double(extractBetween(a{1},"cm**2 = "," a0**2"))*a0^2;
                    end
                    
                    if i == 1 && z == 1 && length(Composition.Z) > 1
                        decs_all = zeros(length(Composition.Z),numel(El.data(:,1)),numel(E0));
                    end
                    if length(Composition.Z) > 1
                        decs_all(z,:,i) = El.data(:,4)*a0^2;
                    end

                    if z == length(Composition.Z) && length(Composition.Z) > 1
                        data.x = El.data(:,1)/180*pi;
                        data.y = sum(decs_all(:,:,i).*Composition.index')/sumweights;
                        data.sigma_el = sum(sigma_el(:,i).*Composition.index')/sumweights;
                        data.sigma_tr1 = sum(sigma_tr1(:,i).*Composition.index')/sumweights;
                        Res(i) = data;
                    elseif length(Composition.Z) == 1
                        data.x = El.data(:,1)/180*pi;
                        if isMolecule
                            data.y = El.data(:,4)*1e16;
                        else
                            data.y = El.data(:,4)*a0^2;
                        end
                        data.sigma_el = sigma_el(1,i);
                        data.sigma_tr1 = sigma_tr1(1,i);
                        Res(i) = data;
                    end                                    
                end
            end

            for i = 5:numel(E0)
                f_name_el=fullfile(dir_Elsepa,...
                    ['dcs_' strrep(strrep(num2str(E0(i), '%1.3e'),'.','p'),'+0','0') '.dat']);
                delete(f_name_el);
            end
        end
    end
    
end

