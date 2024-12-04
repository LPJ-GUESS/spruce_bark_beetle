% Reads the data in the ouput files from logread.m, calibrate all the
% 7^7 = 823543 tested models, calcalte statistics for the models and stors
% the results

clear
close all hidden

% Read textfile with the proportion (weight) of the different stand and age types
weights = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\standweightNy2.txt','\t',1,0);

% Read the list of simulated grid cells  (longitude, latitude)
gridfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\gridlistNy2Clean.txt';
gridlist = dlmread(gridfile,'\t');
gridcells = size(gridlist,1);

% General settings for the simulation
stands = 12;
patches = 5;
ifts = 1;
%years = 2019-1950;
years = 2020-1950;
vars = 13;

% Start and end grid cell number for the County/CountryPart/Country
% included in the simulation. Central Sweden counties (11-18) were included
% in the simulations but not in the analysis
cellregion = [
	1	5	%	1	Skane
	6	9	%	2	Halland
	10	13	%	3	Kronoberg
	14	16	%	4	Blekinge
	17	23	%	5	Kalmar
	24	30	%	6	Alvsborg
	31	36	%	7	Jonkoping
	37	39	%	8	Goteborg
	40	48	%	9	Ostergotland
	49	55	%	10	Skaraborg
	56	60	%	11	Sodermanland
	61	73	%	12	Varmland
	74	78	%	13	Orebro
	79	82	%	14	Stockholm
	83	86	%	15	Vastmanland
	87	91	%	16	Uppsala
	92	111	%	17	Dalarna
	112	124	%	18	Gavleborg
	125	163	%	19	Austria
	164	167	%	20	Alsace
	168	179	%	21	Champagne e Ardenne
	180	187	%	22	Franche Comté
	188	197	%	23	Lorraine
	198	214	%	24	Rhones Alpes
	215	222	%	25	SwissLowland
	223	231	%	26	SwissMountains
    ];
regions = size(cellregion,1);

% The seven parameters for which seven levels were tested in the calibration
params = [
    0.005 0.012 0.025 0.05 0.12 0.25 0.5%k1
    0.15 0.2 0.25 0.35 0.4 0.45 0.5%k2
    -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1%minPIgc
    0.25 0.33333 0.5 1 2 3 4 %previous year weight of wscal
    0.005 0.009 0.017 0.03 0.05 0.08 0.15%k3
    0.25 0.3 0.35 0.4 0.5 0.65 0.8%k4
    5 6 8 10 14 20 27%k5
    ];

% Fixed background level of dead trees available for SBB
base_bm = 0.0001;

% Loading of spruce forest damage fraction (volume) by storm and SBB
load('Marine_obs_damage_filled_SDI.mat'); % obs_damage2ny(region,year,type), 31 years (1990-2020),
% 26 regions (1-10 S Swe, 11-18 C Swe, 19 Aus, 20-24 Fr, 25-26 Swi),2 types (1=storm 2=SBB)
% % For not including Switzerland and Austria data 2011-2019
% obs_damage2ny(:,22:31,:) = NaN;

% Load the Matlab file from the logread script
load('logfileresultsNy2withSSCfix.mat');% !!with/without change 1
%load('logfileresultsNy2withoutSSC.mat');
% inputdata(gridcell,stand,patch,year,var)
% vars: 1=insect_litter, 2=insect_stem, 3=wscal, 4=dam_wood, 5=cmass_wood,
% 6=asp, 7=asp30
% initPI(gridcell,stand,patch)

% Reading of LPJ.GUESS wscal output file
mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\fixedwscalNy2_withSalvSanCut\mwscal_BNE.out','',1,3);% !!with/without change 2
mwscal = reshape(mwscal,[120 gridcells 12]);
mwscal = permute(mwscal,[3 1 2]);
wscalm_current = permute(mean(mwscal(5:7,51:120,:)),[3 2 1]);
wscalm_previous = permute(mean(mwscal(5:7,50:119,:)),[3 2 1]);

% Load previos k0 for reduced number of  calibration iterations
load('IpsCallResNy2withSSC_2010_bak.mat');
k0old = results(:,4);

% Reconstruction the parameter space to two dimensions
paralist = NaN(7^7,7);
startrow = 1;
for k1 = 1:7
    for k2 = 1:7
        for minPIgc = 1:7
            for prevYwscalm = 1:7
                for k3 = 1:7
                    for k4 = 1:7
                        for k5 = 1:7
                            paralist(startrow,1) = params(1,k1);
                            paralist(startrow,2) = params(2,k2);
                            paralist(startrow,3) = params(3,minPIgc);
                            paralist(startrow,4) = params(4,prevYwscalm);% Changed to prevYwscalm from base_bm
                            paralist(startrow,5) = params(5,k3);
                            paralist(startrow,6) = params(6,k4);
                            paralist(startrow,7) = params(7,k5);
                            startrow = startrow + 1;
                        end
                    end
                end
            end
        end
    end
end
paras = size(paralist,1);

max_insect_mort = 0.75;

% Settings for salvage and sanitary cuttting
salvcut = true; % Whether to apply salvage and sanitary cutting !!with/without change 4
sankillfrac = 0.25; % Fraction of bark beetles killed in sanitary counter-measures
sanmin = 0.01; % Min fraction of killed spruce biomass for sanitary cutting to take place

progresscheck = 0;
minlim = 0.8; % Minimum limit of relative deviation for making the final calibration loop
maxlim = 1.3; % Maximum limit of relative deviation for making the final calibration loop
for callsetting = 1:2 %1 includes 2011-2019 calibration data for Austria and Switzerland
    
    results = NaN(paras,4); % 1=R2, 2=RMSE, 3=Bias, 4=k0
    % resreg = NaN(paras,8,3);
    % % 1Norway, 2SSweden, 3CSweden, 4Italy, 5Switzerland, 6Austria, 7France, 8Denmark
    resreg = NaN(paras,4,3);
    % 1SSweden, 2SwitzerlandCSweden, 3Austria, 4France
    % 1=R2, 2=RMSE, 3=Bias
    
    for para = 1:paras

        progresscheck = progresscheck +1;
        if progresscheck > 0.01*paras
            progresscheck = 0;
            callsetting
            progress = para/paras*100
        end

        k0_new = k0old(para);
        withinlimit = false;
        
        while ~withinlimit
        
            %Run the stand-allone SBB modell
            patchPI = initPI;
            mortality_insect = NaN(gridcells,stands,patches,years);
            for year = 2:years
                
                wscalmean = (mean(wscaldata(:,:,:,year,4:6),5) * paralist(para,4) +...
                    mean(wscaldata(:,:,:,year,1:3),5))  / (1 + paralist(para,4));
                droughtR = min(2,(1-wscalmean) * paralist(para,7));

                phenR = inputdata(:,:,:,year,6) * 3.33 ./ (inputdata(:,:,:,year,7) +75);
                patchR = min(1, max((-3.8 - paralist(para,3)), -log(patchPI ./...
                    (inputdata(:,:,:,year,1) + base_bm) * paralist(para,5)) *...
                    paralist(para,6)));
                gcPI = sum(sum(patchPI .* repmat(weights(:,1:12),[1 1 patches]) / patches,3),2);
                gcR = min(1, max(paralist(para,3),-log(gcPI * paralist(para,1)) * paralist(para,2)));
                totR = droughtR + phenR + patchR + repmat(gcR,[1 stands patches]);
                patchPI = exp(totR) .* patchPI;
                % Reduce population if sanitary cutting is applied
                patchPI = patchPI .* (1 - sankillfrac.*...
                    ((patchPI .* k0_new ./ inputdata(:,:,:,year,2)) > sanmin) .*...
                    (inputdata(:,:,:,year,2) > 0) .* salvcut);

                mortality_insect(:,:,:,year) = min(max_insect_mort * inputdata(:,:,:,year,2), k0_new .* patchPI);
            end

            % Summarize for region
            mod_damage = NaN(regions,years-39);%1990-2020
            for region = 1:regions
                startcell = cellregion(region,1);
                endcell = cellregion(region,2);
                cells = endcell - startcell + 1;
                st_weights = weights(startcell:endcell,1:12);
                gc_weights = weights(startcell:endcell,13);
                gc_weights = gc_weights ./ sum(gc_weights);
                all_weights = st_weights .* repmat(gc_weights,[1 12]) ./ patches;
                all_weights = repmat(all_weights,[1 1 patches years-39]);%1990-2020
                reg_ipsmort = mortality_insect(startcell:endcell,:,:,40:end) .* all_weights;%1990-2020
                reg_ipsmort = permute(sum(sum(sum(reg_ipsmort,1))),[1 4 2 3]);
                reg_cmass = inputdata(startcell:endcell,:,:,40:end,5) .* all_weights;%1990-2020
                reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
                mod_damage(region,:) = reg_ipsmort ./ reg_cmass;
            end

            % Adjust k0 based on relative mean deviation from maximum SBB damage
            maxmod = max(mod_damage,[],2);
            maxobs = max(obs_damage2ny(:,:,2),[],2);
            %Include S-Sweden, France, Switzerland and Austria with equal weight
            adjust = (mean(maxobs(1:10)) ./ mean(maxmod(1:10)) + mean(maxobs(20:24)) ./ mean(maxmod(20:24))...
                + mean(maxobs(25:26)) ./ mean(maxmod(25:26)) + maxobs(19) ./ maxmod(19)) / 4;
            k0_new = k0_new * adjust;
            
            % Terminate this loop and make a final run calibration loop if within limit 
            if adjust > minlim && adjust < maxlim
                withinlimit = true;
            end

        end
        
        
        % Final run
        patchPI = initPI;
        mortality_insect = NaN(gridcells,stands,patches,years);
        for year = 2:years

            wscalmean = (mean(wscaldata(:,:,:,year,4:6),5) * paralist(para,4) +...
                mean(wscaldata(:,:,:,year,1:3),5))  / (1 + paralist(para,4));
            droughtR = min(2,(1-wscalmean) * paralist(para,7));

            phenR = inputdata(:,:,:,year,6) * 3.33 ./ (inputdata(:,:,:,year,7) +75);
            patchR = min(1, max((-3.8 - paralist(para,3)), -log(patchPI ./...
                (inputdata(:,:,:,year,1) + base_bm) * paralist(para,5)) *...
                paralist(para,6)));
            gcPI = sum(sum(patchPI .* repmat(weights(:,1:12),[1 1 patches]) / patches,3),2);
            gcR = min(1, max(paralist(para,3),-log(gcPI * paralist(para,1)) * paralist(para,2)));
            totR = droughtR + phenR + patchR + repmat(gcR,[1 stands patches]);
            patchPI = exp(totR) .* patchPI;
            % Reduce population if sanitary cutting is applied
            patchPI = patchPI .* (1 - sankillfrac.*...
                ((patchPI .* k0_new ./ inputdata(:,:,:,year,2)) > sanmin) .*...
                (inputdata(:,:,:,year,2) > 0) .* salvcut);

            mortality_insect(:,:,:,year) = min(max_insect_mort * inputdata(:,:,:,year,2), k0_new .* patchPI);
        end

        % Summarize for region
        mod_damage = NaN(regions,years-39);%1990-2020
        r2region = NaN(regions,1);
        for region = 1:regions
            startcell = cellregion(region,1);
            endcell = cellregion(region,2);
            cells = endcell - startcell + 1;
            st_weights = weights(startcell:endcell,1:12);
            gc_weights = weights(startcell:endcell,13);
            gc_weights = gc_weights ./ sum(gc_weights);
            all_weights = st_weights .* repmat(gc_weights,[1 12]) ./ patches;
            all_weights = repmat(all_weights,[1 1 patches years-39]);%1990-2020
            reg_ipsmort = mortality_insect(startcell:endcell,:,:,40:end) .* all_weights;%1990-2020
            reg_ipsmort = permute(sum(sum(sum(reg_ipsmort,1))),[1 4 2 3]);
            reg_cmass = inputdata(startcell:endcell,:,:,40:end,5) .* all_weights;%1990-2020
            reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
            mod_damage(region,:) = reg_ipsmort ./ reg_cmass;

            regcorr = corrcoef((mod_damage(region,:))',(obs_damage2ny(region,:,2))','Rows','complete');
            r2region(region) = regcorr(1,2) ^2;
        end
        
        % Calculate statistics, R2, RMSE, Bias
        sqrdiffs = (mod_damage - obs_damage2ny(:,:,2)) .^2;
        rmse = mean(sqrdiffs,2,'omitnan') .^0.5;
        bias = mean(mod_damage - obs_damage2ny(:,:,2),2,'omitnan');

        % Summarize for large regions countries
        % 1=R2, 2=RMSE, 3=Bias
        resreg(para,1,1) = mean(r2region(1:10));%Southern Sweden
        resreg(para,2,1) = mean(r2region(25:26));%Switzerland
        resreg(para,3,1) = r2region(19);%Austria
        resreg(para,4,1) = mean(r2region(20:24));%France

        resreg(para,1,2) = mean(rmse(1:10));
        resreg(para,2,2) = mean(rmse(25:26));
        resreg(para,3,2) = rmse(19);
        resreg(para,4,2) = mean(rmse(20:24));

        resreg(para,1,3) = mean(bias(1:10));
        resreg(para,2,3) = mean(bias(25:26));
        resreg(para,3,3) = bias(19);
        resreg(para,4,3) = mean(bias(20:24));

        results(para,4) = k0_new; % 1=R2, 2=RMSE, 3=Bias, 4=k0
    end

    % Average over S Sweden and NE France, Switzerland and Austria
    results(:,1) = (resreg(:,1,1) + resreg(:,4,1) + resreg(:,3,1) + resreg(:,2,1)) / 4;% R2 not weighted
    invsum = 1 / mean(maxobs(1:10)) + 1 / mean(maxobs(20:24)) + 1 / mean(maxobs(25:26)) + 1 / maxobs(19);
    results(:,2) = (resreg(:,1,2) ./ mean(maxobs(1:10)) +...
        resreg(:,2,2) ./ mean(maxobs(25:26)) +...
        resreg(:,4,2) ./ mean(maxobs(20:24)) +...
        resreg(:,3,2) ./ maxobs(19)) / invsum;
    results(:,3) = (resreg(:,1,3) ./ mean(maxobs(1:10)) +...
        resreg(:,2,3) ./ mean(maxobs(25:26)) +...
        resreg(:,4,3) ./ mean(maxobs(20:24)) +...
        resreg(:,3,3) ./ maxobs(19)) / invsum;

    if callsetting == 1
        save('IpsCallResNy2withSSCfix_2019.mat','results','resreg','paralist');% !!with/without change 5
        
        % For not including Switzerland and Austria data 2011-2019 in the second run
        obs_damage2ny(:,22:31,:) = NaN;
    else
        save('IpsCallResNy2withSSCfix_2010.mat','results','resreg','paralist');% !!with/without change 6
    end
end

