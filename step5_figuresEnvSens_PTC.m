% Makes figures of stand-allone vs LPJ-GUESS run with the "best" calibrated
% model, and figure of Predsiposing, Triggering and Contributing factors

clear
close all hidden

    
params = [
    0.005 0.012 0.025 0.05 0.12 0.25 0.5%k1
    0.15 0.2 0.25 0.35 0.4 0.45 0.5%k2
    -2.5 -2.25 -2 -1.75 -1.5 -1.25 -1%minPIgc
    0.25 0.33333 0.5 1 2 3 4 %previous year weight of wscal
    0.005 0.009 0.017 0.03 0.05 0.08 0.15%k3
    0.25 0.3 0.35 0.4 0.5 0.65 0.8%k4
    5 6 8 10 14 20 27%k5
    ];

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

stands = 12;
patches = 5;
years = 2020-1950;

max_insect_mort = 0.75;
base_bm = 0.0001; % Still use a very low number to avoid devision by 0

salvcut = true;
sankillfrac = 0.25; % Fraction of bark beetles killed in sanitary counter-measures
sanmin = 0.01; % Min fraction of killed spruce biomass for sanitary cutting to take place

weights = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\standweightNy2.txt','\t',1,0);

load('Marine_obs_damage_filled_SDI.mat'); % obs_damage2ny(region,year,type), 31 years (1990-2020),
% 26 regions (1-10 S Swe, 11-18 C Swe, 19 Aus, 20-24 Fr, 25-26 Swi),2 types (1=storm 2=Ips)
obs_damage = NaN(26,60,2);
obs_damage(:,30:60,:) = obs_damage2ny;

timedata = NaN(regions,years-10,7); % ,reg,year,serie
% (1=storm_obs,2=storm_call,3=storm_call2+,4=insect_obs,5=insect_call,6=insect_call2+,7=defaultcallibrated
timedata(:,:,1) = obs_damage(:,:,1);
timedata(:,:,4) = obs_damage(:,:,2);
timeplotdata = NaN(4,years-10,7); % ,large_reg,year,serie

timedata2 = NaN(regions,years-10,9); % ,reg,year,serie
% 1=cmass call, 2=cmass call +2, 3=cmass default, 4=insect_stem call, 5=insect_stem +2,
% 6=insect_stem default, 7=abs damage call, 8=abs damage call +2, 9=abs damage default
timeplotdata2 = NaN(4,years-10,9); % ,large_reg,year,serie

timedata3 = NaN(regions,years-10,10); % ,reg,year,serie
% 1=droughtR call, 2=gcR call, 3=patchR call, 4=droughtR def, 5=gcR def, 6=patchR def
% 7=patchPI call, 8=patchPI def, 9=totR call, 10=totR def
timeplotdata3 = NaN(4,years-10,10); % ,large_reg,year,serie

%Modelled result of callibrated model
load('logfileresultsNy2callibrated.mat');
%load('logfileresultsNy2withSSCfix.mat');
load('guessresultsNy2callibrated.mat');

damfrac_mod(isnan(damfrac_mod)) = 0;

mortality_storm = damfrac_mod(:,:,:,:,1);
mortality_insect = damfrac_mod(:,:,:,:,2);
insect_stem = damfrac_mod(:,:,:,:,3);

% Summarize for region
for region = 1:regions
    startcell = cellregion(region,1);
    endcell = cellregion(region,2);
    cells = endcell - startcell + 1;
    st_weights = weights(startcell:endcell,1:12);
    gc_weights = weights(startcell:endcell,13);
    gc_weights = gc_weights ./ sum(gc_weights);
    all_weights = st_weights .* repmat(gc_weights,[1 12]) ./ patches;
    all_weights = repmat(all_weights,[1 1 patches years-10]);
    reg_ipsmort = mortality_insect(startcell:endcell,:,:,11:end) .* ...
        insect_stem(startcell:endcell,:,:,11:end) .*all_weights;
    reg_ipsmort = permute(sum(sum(sum(reg_ipsmort,1))),[1 4 2 3]);
    reg_stormmort = mortality_storm(startcell:endcell,:,:,11:end) .* ...
        inputdata(startcell:endcell,:,:,11:end,5) .*all_weights;
    reg_stormmort = permute(sum(sum(sum(reg_stormmort,1))),[1 4 2 3]);
    reg_cmass = inputdata(startcell:endcell,:,:,11:end,5) .* all_weights;
    reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
    reg_insect_stem = insect_stem(startcell:endcell,:,:,11:end) .*all_weights;
    reg_insect_stem = permute(sum(sum(sum(reg_insect_stem,1))),[1 4 2 3]);

    timedata(region,:,2) = reg_stormmort ./ reg_cmass;
    timedata(region,:,5) = reg_ipsmort ./ reg_cmass;
    
    timedata2(region,:,1) = reg_cmass;
    timedata2(region,:,4) = reg_insect_stem;
    timedata2(region,:,7) = reg_ipsmort;
    kolla = 1;
end
    
%Modelled result of callibrated model with +2 degrees
load('logfileresultsNy2callibrated2p.mat');
%load('logfileresultsNy2withoutSSCfix.mat');


damfrac_mod(isnan(damfrac_mod)) = 0;

mortality_storm = damfrac_mod(:,:,:,:,1);
mortality_insect = damfrac_mod(:,:,:,:,2);
insect_stem = damfrac_mod(:,:,:,:,3);

% Summarize for region
for region = 1:regions
    startcell = cellregion(region,1);
    endcell = cellregion(region,2);
    cells = endcell - startcell + 1;
    st_weights = weights(startcell:endcell,1:12);
    gc_weights = weights(startcell:endcell,13);
    gc_weights = gc_weights ./ sum(gc_weights);
    all_weights = st_weights .* repmat(gc_weights,[1 12]) ./ patches;
    all_weights = repmat(all_weights,[1 1 patches years-10]);
    reg_ipsmort = mortality_insect(startcell:endcell,:,:,11:end) .* ...
        insect_stem(startcell:endcell,:,:,11:end) .*all_weights;
    reg_ipsmort = permute(sum(sum(sum(reg_ipsmort,1))),[1 4 2 3]);
    reg_stormmort = mortality_storm(startcell:endcell,:,:,11:end) .* ...
        inputdata(startcell:endcell,:,:,11:end,5) .*all_weights;
    reg_stormmort = permute(sum(sum(sum(reg_stormmort,1))),[1 4 2 3]);
    reg_cmass = inputdata(startcell:endcell,:,:,11:end,5) .* all_weights;
    reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
    reg_insect_stem = insect_stem(startcell:endcell,:,:,11:end) .*all_weights;
    reg_insect_stem = permute(sum(sum(sum(reg_insect_stem,1))),[1 4 2 3]);

    timedata(region,:,3) = reg_stormmort ./ reg_cmass;
    timedata(region,:,6) = reg_ipsmort ./ reg_cmass;
    
    timedata2(region,:,2) = reg_cmass;
    timedata2(region,:,5) = reg_insect_stem;
    timedata2(region,:,8) = reg_ipsmort;
end


% Parameter-combination number for the "best" model  calibrated without Aus
% and Swi 2011-2019 data with SSC, run with SSC
para = 108382;
%para = 393249; % Default parameter set
k0_new = 0.005124114; % New wscalfix and initfix callibrated value

load('logfileresultsNy2withSSCfix.mat');
gridcells = size(inputdata,1);

modelstate = NaN(gridcells,stands,patches,years,7);
% 1=droughtR,2=patchR,3=gcPI,4=gcR,5=patchR,6=patchPI,7=mortality_insect

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
    patchPI = patchPI .* (1 - sankillfrac.*...
        ((patchPI .* k0_new ./ inputdata(:,:,:,year,2)) > sanmin) .*...
        (inputdata(:,:,:,year,2) > 0) .* salvcut);

    mortality_insect(:,:,:,year) = min(max_insect_mort * inputdata(:,:,:,year,2), k0_new .* patchPI);

    modelstate(:,:,:,year,1) = droughtR;
    modelstate(:,:,:,year,2) = phenR;
    modelstate(:,:,:,year,3) = repmat(gcPI,[1 stands patches]);
    modelstate(:,:,:,year,4) = repmat(gcR,[1 stands patches]);
    modelstate(:,:,:,year,5) = patchR;
    modelstate(:,:,:,year,6) = patchPI;
    modelstate(:,:,:,year,7) = mortality_insect(:,:,:,year);
end

insect_stem = damfrac_mod(:,:,:,:,3);

% Summarize for region
mod_damage = NaN(regions,years-10);
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
    reg_cmass = inputdata(startcell:endcell,:,:,11:end,5) .* all_weights;
    reg_cmass = permute(sum(sum(sum(reg_cmass,1))),[1 4 2 3]);
    mod_damage(region,:) = reg_ipsmort ./ reg_cmass;
    reg_insect_stem = insect_stem(startcell:endcell,:,:,11:end) .*all_weights;
    reg_insect_stem = permute(sum(sum(sum(reg_insect_stem,1))),[1 4 2 3]);

    timedata(region,:,7) = mod_damage(region,:);
    
    timedata2(region,:,3) = reg_cmass;
    timedata2(region,:,6) = reg_insect_stem;
    timedata2(region,:,9) = reg_ipsmort;
    
    tempR = guessresults(startcell:endcell,:,:,11:end,1) .*all_weights; %droughtR callibrated
    timedata3(region,:,1) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = guessresults(startcell:endcell,:,:,11:end,4) .*all_weights; %gcR callibrated
    timedata3(region,:,2) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = guessresults(startcell:endcell,:,:,11:end,5) .*all_weights; %patchR callibrated
    timedata3(region,:,3) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = modelstate(startcell:endcell,:,:,11:end,1) .*all_weights; %droughtR default
    timedata3(region,:,4) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = modelstate(startcell:endcell,:,:,11:end,4) .*all_weights; %gcR default
    timedata3(region,:,5) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = modelstate(startcell:endcell,:,:,11:end,5) .*all_weights; %patchR default
    timedata3(region,:,6) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = guessresults(startcell:endcell,:,:,11:end,6) .*all_weights; %patchPI callibrated
    timedata3(region,:,7) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    tempR = modelstate(startcell:endcell,:,:,11:end,6) .*all_weights; %patchPI default
    timedata3(region,:,8) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    
    tempR = guessresults(startcell:endcell,:,:,11:end,2) .*all_weights; %totalR callibrated
    timedata3(region,:,9) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    timedata3(region,:,9) = timedata3(region,:,9) + timedata3(region,:,1) +...
        timedata3(region,:,2) + timedata3(region,:,3);
    tempR = modelstate(startcell:endcell,:,:,11:end,2) .*all_weights; %totalR default
    timedata3(region,:,10) = permute(sum(sum(sum(tempR,1))),[1 4 2 3]);
    timedata3(region,:,10) = timedata3(region,:,10) + timedata3(region,:,4) +...
        timedata3(region,:,5) + timedata3(region,:,6);
end


% Summarize for large regions/countries
timeplotdata(1,:,:) = mean(timedata(1:10,:,:));%Southern Sweden
timeplotdata(2,:,:) = mean(timedata(25:26,:,:));%Switzerland
timeplotdata(3,:,:) = timedata(19,:,:);%Austria
timeplotdata(4,:,:) = mean(timedata(20:24,:,:));%France

timeplotdata2(1,:,:) = mean(timedata2(1:10,:,:));%Southern Sweden
timeplotdata2(2,:,:) = mean(timedata2(25:26,:,:));%Switzerland
timeplotdata2(3,:,:) = timedata2(19,:,:);%Austria
timeplotdata2(4,:,:) = mean(timedata2(20:24,:,:));%France

timeplotdata3(1,:,:) = mean(timedata3(1:10,:,:));%Southern Sweden
timeplotdata3(2,:,:) = mean(timedata3(25:26,:,:));%Switzerland
timeplotdata3(3,:,:) = timedata3(19,:,:);%Austria
timeplotdata3(4,:,:) = mean(timedata3(20:24,:,:));%France

colors = [0 0 0%Black
    0.1 0.25 1%Dark blue
    0.6 0.3 0.07%Dark brown
    0.5 0.5 0.5%Grey
    0.3 0.80 1%Bright blue
    0.93 0.65 0.43%Light brown
    1 0.65 0.8%Pink
    ];

styles = {'-'
    '-'
    '--'
    '--'
    '-'
    '-'
    '--'
    '--'
    '-'};

markers = 'sdosdo^';

xstart = 0.06;
ystart = 0.09;
xwidth = 0.88;
yheight = 0.88;

xdata = 1961:2020;

legtext = {
    'Storm obs'
    'Storm sim'
    'Storm sim+2'
    'SBB obs'
    'SBB sim'
    'SBB sim+2'
    'SBB def'};


path = 'C:\Users\Fredrik.Lagergren\Documents\BioticMod\manus\Testing\DefaultTimeplots_';
regtext = {
    'S_Swe'
    'Switz'
    'Austria'
    'France'};


%Timeplots
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
    hold(axel1,'all');
    xlim(axel1,[1990 2020]);
    
    yyaxis left; 
    ylabel(axel1,'% damaged storm');
    hold(axel1,'all');
    for serie = 1:3
        ydata = permute(timeplotdata(large_reg,:,serie),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(serie,:),'LineWidth',2,'LineStyle','-',...
            'Marker',markers(serie),'MarkerEdgeColor','none','MarkerFaceColor',colors(serie,:));
        hold(axel1,'all');
    end
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    hold(axel1,'all');
    
    yyaxis right; 
    ylabel(axel1,'% damaged SBB');
    hold(axel1,'all');
    for serie = 4:7
        ydata = permute(timeplotdata(large_reg,:,serie),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(serie,:),'LineWidth',3,'LineStyle','-',...
            'Marker',markers(serie),'MarkerEdgeColor','none','MarkerFaceColor',colors(serie,:));
        hold(axel1,'all');
    end
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    hold(axel1,'all');

    legend(axel1,legtext,...
        'FontSize',16,'Location','northwest','NumColumns',2);
    thefile = sprintf('%s',path,'Timecomp_wscalfix',regtext{large_reg},'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end


% Timeplots predesposing, triggering and contributing factors - calibrated data
timeplotdata4 = NaN(4,years-10,4);
timeplotdata4(:,:,1) = timeplotdata(:,:,1) * 100; % Relative storm damage
timeplotdata4(:,:,1) = timeplotdata4(:,:,1) / 5; % To have common scale with biomass
%timeplotdata4(:,:,1) = timeplotdata2(:,:,9); % Absolute storm damage
timeplotdata4(:,:,2) = timeplotdata2(:,:,4);
timeplotdata4(:,:,3) = timeplotdata3(:,:,10) - timeplotdata3(:,:,5) - timeplotdata3(:,:,6);
timeplotdata4(:,3:60,4) = (timeplotdata4(:,1:58,3) + timeplotdata4(:,2:59,3) + timeplotdata4(:,3:60,3)) / 3;
legtext = {
    'Storm damage'
    'Spruce > 15 cm'
    'f(phen.) + f(water stress)'
    'f(ph.) + f(w. s.) 3-yr mean'
    };

colors = [
    0.1 0.25 1%Dark blue
    1 0.65 0.8%Pink
    0.6 0.3 0.07%Dark brown
    0.93 0.65 0.43%Bright brown
    ];
markers = 'sdoo^';

xstart = 0.12;
xwidth = 0.81;
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
    hold(axel1,'all');
    xlim(axel1,[1990 2020]);
    
    yyaxis left; 
    ylabel(axel1,'Biomass (kg C m^-^2)');
    hold(axel1,'all');
    for serie = 1:2
        ydata = permute(timeplotdata4(large_reg,:,serie),[2 1]);
        plot(axel1,xdata, ydata,'Color',colors(serie,:),'LineWidth',3,'LineStyle','-',...
            'Marker',markers(serie),'MarkerEdgeColor','none','MarkerFaceColor',colors(serie,:));
        hold(axel1,'all');
    end
    ylim(axel1,[0 5]);
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    hold(axel1,'all');
    
    yyaxis right; 
    ylabel(axel1,'R component sum');
    hold(axel1,'all');
    for serie = 3:4
        ydata = permute(timeplotdata4(large_reg,:,serie),[2 1]);
        if serie == 3
            [p,S] = polyfit(xdata(30:60)',ydata(30:60),1);
        end
        plot(axel1,xdata, ydata,'Color',colors(serie,:),'LineWidth',3,'LineStyle','-',...
            'Marker',markers(serie),'MarkerEdgeColor','none','MarkerFaceColor',colors(serie,:));
        hold(axel1,'all');
    end
    plot(axel1,[1990 2020], [1990*p(1)+p(2) 2020*p(1)+p(2)],'Color',colors(3,:),...
        'LineWidth',2,'LineStyle','--','Marker','none');
    hold(axel1,'all');
    ylim(axel1,[0 3]);
    yticks(axel1,0:0.6:3);
    hold(axel1,'all');
    ax = gca;
    ax.YColor = 'k';
    hold(axel1,'all');

    legend(axel1,legtext,...
        'FontSize',16,'Location','northwest','NumColumns',2);
    axel2 = axes('Position',[xstart-0.06 ystart 0 yheight],'Box','off','XGrid','off','YGrid','off',...
        'FontSize',16,'Parent',figure1);
    ylabel(axel2,'Damaged storm (%)');
    ylim(axel2,[0 25]);
   
    thefile = sprintf('%s',path,'Sensitivity',regtext{large_reg},'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end

