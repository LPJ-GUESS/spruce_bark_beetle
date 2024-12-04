% Based on sorting of the stats, the best and 50 best models are found for
% each region/country and combined. Results are written to tables exported
% to Excel and then further to Word for presentation. Time plots showing
% observation, default simulation results and with the "best" combined
% stats for the four region/countries are saved as emf-files

clear
close all hidden

% Load results from callibration.m
load('IpsCallResNy2withSSCfix_2010.mat'); % !!with/without + 2010/2019 change 1
% results(para,resvar) % 1=R2, 2=RMSE, 3=Bias, 4=k0)
% resreg(para,reg,resvar)
% paralist(para,parvar) 1=k1, 2=k2, 3=minPIgc, 4=base_bm, 5=k3, 6=k4, 7=k5

path = 'C:\Users\Fredrik.Lagergren\Documents\BioticMod\IPSEuropeSens\WithSSC\fixwscal\'; % !!with/without change 2

regtext = {
    'All4'
    'S_Swe'
    'Switz'
    'Austria'
    'France'};

stattext = {
    'R2_'
    'RMSE_'
    'Bias_'
    'Combind'};

bestpars = NaN(5,9*4);
table_stats = NaN(5,9);


%Find lowest normalized sum of 1-R2, rmse and |bias|
%Deviation from mean divided by range
restotal = (1 - results(:,1)) - (1 - mean(results(:,1))) / (max(1 - results(:,1)) - min(1 - results(:,1))) +...%R2
    (results(:,2) - mean(results(:,2))) / (max(results(:,2)) - min(results(:,2))) +...%RMSE
    (abs(results(:,3)) - mean(abs(results(:,3)))) / (max(abs(results(:,3))) - min(abs(results(:,3))));%Bias
[~,index] = sortrows(restotal);
ind4select = index(1);
ind50best = index(1:50);
ind_reg = NaN(4,1); % Number of the best combined parameterization
ind_reg50best = NaN(50,4);
for reg = 1:4
    restotal = (1 - resreg(:,reg,1)) - (1 - mean(resreg(:,reg,1))) / (max(1 - resreg(:,reg,1)) - min(1 - resreg(:,reg,1))) +...%R1
        (resreg(:,reg,2) - mean(resreg(:,reg,2))) / (max(resreg(:,reg,2)) - min(resreg(:,reg,2))) +...%RMSE
        (abs(resreg(:,reg,3)) - mean(abs(resreg(:,reg,3)))) / (max(abs(resreg(:,reg,3))) - min(abs(resreg(:,reg,3))));%Bias
    [~,index] = sortrows(restotal);
    ind_reg(reg) = index(1);
    ind_reg50best(:,reg) = index(1:50);
end
    
    
ind_def = 393249; % Number of the parameterization used for the undelying Guess simulations (4 3 3 6 4 4 3) Ny2SSC


for stattype = 1:4% 1 = R2, 2 = RMSE, 3 = Bias, 4 = combined

    for reg = 0:4

        if stattype < 4
            if reg == 0
                resdata = results(:,stattype);

            else
                resdata = resreg(:,reg,stattype);
            end
        else
            if reg == 0% All regions combined
                %resdata = indexsum4sel;
                resdata = results(:,3); %This will be the same as bias
            else
                %resdata = indexsumreg(:,reg);
                resdata = resreg(:,reg,3); %This will be the same as bias
            end
        end

        if stattype == 1% R2, highest value best
            [~,index] = sortrows(resdata,'descend');
        else
            resdata = resdata * 100;
            if stattype == 3% Bias, lowest absolute value best
                [~,index] = sortrows(abs(resdata));
            else% 2RMSE 4Combined, lowest value best
                [~,index] = sortrows(resdata);
            end
        end
        index = index(1);
        if stattype == 1 && reg == 0
            %ind4select = index;
            table_stats(1,1:3) = results(ind_def,1:3);
            table_stats(1,4:6) = results(ind4select,1:3);
            table_stats(2:5,1:3) = permute(resreg(ind_def,:,:),[2 3 1]);
            table_stats(2:5,4:6) = permute(resreg(ind4select,:,:),[2 3 1]);
        end
        coladd = (stattype-1) * 9;
        if stattype < 4
            bestpars(reg+1,1+coladd:7+coladd) = paralist(index,:); % parameters
            bestpars(reg+1,8+coladd) = results(index,4); % k0
            bestpars(reg+1,9+coladd) = resdata(index); % stat r2
        else
            if reg == 0
                bestpars(reg+1,1+coladd:7+coladd) = paralist(ind4select,:); % parameters
                bestpars(reg+1,8+coladd) = results(ind4select,4); % k0
                bestpars(reg+1,9+coladd) = resdata(ind4select); % stat r2
               
            else
                bestpars(reg+1,1+coladd:7+coladd) = paralist(ind_reg(reg),:); % parameters
                bestpars(reg+1,8+coladd) = results(ind_reg(reg),4); % k0
                bestpars(reg+1,9+coladd) = resdata(ind_reg(reg)); % stat
            end
        end
            
        if stattype == 1 && reg > 0
            table_stats(reg+1,7:9) = permute(resreg(ind_reg(reg),reg,:),[2 3 1]);
        end

    end
end


% Plot time graphs

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


load('Marine_obs_damage_filled_SDI.mat'); % obs_damage2ny(region,year,type), 31 years (1990-2020),
% 26 regions (1-10 S Swe, 11-18 C Swe, 19 Aus, 20-24 Fr, 25-26 Swi),2 types (1=storm 2=Ips)

% Extend the serie to 1961-2020 with NaN to have the same period as the
% results of the stand-alone module
obs_damage = NaN(26,60,2);
obs_damage(:,30:60,:) = obs_damage2ny;

load('logfileresultsNy2withSSCfix.mat'); % !!with/without change 4
% inputdata(gridcell,stand,patch,year,var)
% vars: 1=insect_litter, 2=insect_stem, 3=wscal, 4=dam_wood, 5=cmass_wood,
% 6=asp, 7=asp30
% initPI(gridcell,stand,patch)
gridcells = size(inputdata,1);
stands = 12;
patches = 5;
years = 2020-1950;

paras = NaN(6,1);
paras(1) = ind_def;
paras(2) = ind4select;
paras(3:6) = ind_reg;

max_insect_mort = 0.75;

% Fixed background level of dead trees available for SBB
base_bm = 0.0001;

salvcut = true; % Whether to apply salvage and sanitary cutting  % !!with/without change 5
sankillfrac = 0.25; % Fraction of bark beetles killed in sanitary counter-measures
sanmin = 0.01; % Min fraction of killed spruce biomass for sanitary cutting to take place

weights = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\standweightNy2.txt','\t',1,0);

mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\fixedwscalNy2_withSalvSanCut\mwscal_BNE.out','',1,3); % !!with/without change 6
mwscal = reshape(mwscal,[120 gridcells 12]);
mwscal = permute(mwscal,[3 1 2]);
wscalm_current = permute(mean(mwscal(5:7,51:120,:)),[3 2 1]);
wscalm_previous = permute(mean(mwscal(5:7,50:119,:)),[3 2 1]);

timedata = NaN(regions,years-10,4); % ,reg,year,type (1=obs, 2=default, 3=4sel, 4=reg; 5-8 same for storm)
timedata(:,:,1) = obs_damage(:,:,2);
timedata(:,:,5) = obs_damage(:,:,1);

% Run the stand-allone module for the default (1) and "best" models selected
% for all four regions/countries (2) of by specific region/country (3-6)
for combo = 1:6
    para = paras(combo);
    k0_new = results(para,4);

    patchPI = initPI;
    mortality_insect = NaN(gridcells,stands,patches,years);
    mortality_storm = NaN(gridcells,stands,patches,years);
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
        patchPI = patchPI .* (1 - sankillfrac.*...
            ((patchPI .* k0_new ./ inputdata(:,:,:,year,2)) > sanmin) .*...
            (inputdata(:,:,:,year,2) > 0) .* salvcut);

        mortality_insect(:,:,:,year) = min(max_insect_mort * inputdata(:,:,:,year,2), k0_new .* patchPI);
        mortality_storm(:,:,:,year) = inputdata(:,:,:,year,4);
    end

    % Summarize for region
    mod_damage = NaN(regions,years-10);
    mod_stormdamage = NaN(regions,years-10);
    r2region = NaN(regions,1);
    for region = 1:regions
        startcell = cellregion(region,1);
        endcell = cellregion(region,2);
        cells = endcell - startcell + 1;
        st_weights = weights(startcell:endcell,1:12);
        gc_weights = weights(startcell:endcell,13);
        gc_weights = gc_weights ./ sum(gc_weights);
        all_weights = st_weights .* repmat(gc_weights,[1 12]) ./ patches;
        all_weights = repmat(all_weights,[1 1 patches years-10]);
        reg_ipsmort = mortality_insect(startcell:endcell,:,:,11:end) .* all_weights;
        reg_ipsmort = permute(sum(sum(sum(reg_ipsmort,1))),[1 4 2 3]);
        reg_stormmort = mortality_storm(startcell:endcell,:,:,11:end) .* all_weights;
        reg_stormmort = permute(sum(sum(sum(reg_stormmort,1))),[1 4 2 3]);
        reg_cmass = inputdata(startcell:endcell,:,:,11:end,5) .* all_weights;
        reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
        mod_damage(region,:) = reg_ipsmort ./ reg_cmass;
        mod_stormdamage(region,:) = reg_stormmort ./ reg_cmass;
        
        large_reg = NaN;
        switch(region)
           case {1,2,3,4,5,6,7,8,9,10}
              large_reg = 1;%Southern Sweden
           case {25,26}
              large_reg = 2;%Switzerland
           case 19
              large_reg = 3;%Austria
           case {20,21,22,23,24}
              large_reg = 4;%France
        end   
        if combo-2 == large_reg
            timedata(region,:,4) = mod_damage(region,:);
            timedata(region,:,8) = mod_stormdamage(region,:);
        end
        if combo < 3
            timedata(region,:,combo+1) = mod_damage(region,:);
            timedata(region,:,combo+1+4) = mod_stormdamage(region,:);
        end
        
        regcorr = corrcoef((mod_damage(region,:))',(obs_damage(region,:,2))','Rows','complete');
        r2region(region) = regcorr(1,2) ^2;
    end
    
end

% Summarize for large regions/countries
timeplotdata = NaN(4,years-10,8);
timeplotdata(1,:,:) = mean(timedata(1:10,:,:));%Southern Sweden
timeplotdata(2,:,:) = mean(timedata(25:26,:,:));%Switzerland
timeplotdata(3,:,:) = timedata(19,:,:);%Austria
timeplotdata(4,:,:) = mean(timedata(20:24,:,:));%France

stddata = NaN(4,years-10,6); %1-2 +/- obs std ips, 3-4 storm, 5 std ips, 6 std storm
stddata(1,:,1) = timeplotdata(1,:,1) - std(timedata(1:10,:,1));
stddata(1,:,2) = timeplotdata(1,:,1) + std(timedata(1:10,:,1));
stddata(2,:,1) = timeplotdata(2,:,1) - std(timedata(25:26,:,1));
stddata(2,:,2) = timeplotdata(2,:,1) + std(timedata(25:26,:,1));
stddata(4,:,1) = timeplotdata(4,:,1) - std(timedata(20:24,:,1));
stddata(4,:,2) = timeplotdata(4,:,1) + std(timedata(20:24,:,1));
stddata(1,:,3) = timeplotdata(1,:,5) - std(timedata(1:10,:,5));
stddata(1,:,4) = timeplotdata(1,:,5) + std(timedata(1:10,:,5));
stddata(2,:,3) = timeplotdata(2,:,5) - std(timedata(25:26,:,5));
stddata(2,:,4) = timeplotdata(2,:,5) + std(timedata(25:26,:,5));
stddata(4,:,3) = timeplotdata(4,:,5) - std(timedata(20:24,:,5));
stddata(4,:,4) = timeplotdata(4,:,5) + std(timedata(20:24,:,5));

stddata(1,:,5) = std(timedata(1:10,:,1));
stddata(2,:,5) = std(timedata(25:26,:,1));
stddata(4,:,5) = std(timedata(20:24,:,1));
stddata(1,:,6) = std(timedata(1:10,:,5));
stddata(2,:,6) = std(timedata(25:26,:,5));
stddata(4,:,6) = std(timedata(20:24,:,5));

stddata(stddata < 0) = 0;

xdata = 1961:2020;

colors = [0 0 0
    0.1 0.25 1
    0.1 0.3 0.1
    0.8 0.4 0.8
    0.6 0 0
    0.9 0.5 0.1];

markers = 'osd^+x';

xstart = 0.06;
ystart = 0.09;
xwidth = 0.88;
yheight = 0.88;

% Timeplots
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','off','YGrid','off',...
        'FontSize',16,'Parent',figure1);
    %hold(axel1,'all');
    xlim(axel1,[1990 2020]);
    
    yyaxis left; 
    ylabel(axel1,'% damaged storm');
    hold(axel1,'all');

    ydata = permute(timeplotdata(large_reg,:,5),[2 1]) * 100;
    plotdata = NaN(60,2);
    plotdata(:,1) = xdata;
    plotdata(:,2) = ydata;
    plotdata = sortrows(plotdata,2);
    finrows = sum(isfinite(ydata));
    plotdata = plotdata(1:finrows,:);
    plotdata = sortrows(plotdata,1);
    bar(axel1,plotdata(:,1) - 0.10,plotdata(:,2),'BarWidth',0.15,'FaceColor',colors(5,:),'LineStyle','-');
    hold(axel1,'all');
    
    ydata = permute(timeplotdata(large_reg,:,6),[2 1]) * 100;
    bar(axel1,xdata + 0.10,ydata,'BarWidth',0.15,'FaceColor',colors(6,:),'LineStyle','-');
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    hold(axel1,'all');
    
    yyaxis right;
    ylabel(axel1,'% damaged SBB');
    hold(axel1,'all');
    
    for partype = 1:4
        ydata = permute(timeplotdata(large_reg,:,partype),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(partype,:),'LineWidth',2,'LineStyle','-',...
            'Marker',markers(partype),'MarkerEdgeColor','none','MarkerFaceColor',colors(partype,:));
    end
    hold(axel1,'all');
    
    yyaxis right;
    for stdshow = 1:2
        ydata = permute(stddata(large_reg,:,stdshow),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(1,:),'LineWidth',2,'LineStyle',':','Marker','none');
    end
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    
    legend(axel1,{'Storm obs','Storm sim','SBB obs','Default','All4 opt','Region opt'},...
        'FontSize',16,'Location','northwest');
    hold(axel1,'all');
        
    
    thefile = sprintf('%s',path,'Timeplot',regtext{large_reg+1},'_10.emf'); % !!2010/2019 change 7
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end

par50best = paralist(ind50best,:);
min50best = min(par50best);
max50best = max(par50best);

