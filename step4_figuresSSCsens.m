% Runs the stand-alone model and makes figures for the sensitivity of
% including Salvage and Sanetary Cutting (SSC) by running the "best" models
% calibrated with/without SSC with/without SSC setting. Additionaly,
% ackumulated figures are made.

clear
close all hidden

% Headers for the different settings
legtext = {
    'cW10-rW'% Cal. without Aus and Swi 2011-2019 data with SSC, run with SSC
    'cW19-rW'% Cal. with Aus and Swi 2011-2019 data with SSC, run with SSC
    'cW10-rWO'% Cal. without Aus and Swi 2011-2019 data with SSC, run without SSC
    'cW19-rWO'% Cal. with Aus and Swi 2011-2019 data with SSC, run without SSC
    'cWO10-rWO'% Cal. without Aus and Swi 2011-2019 data without SSC, run without SSC
    'cWO19-rWO'% Cal. with Aus and Swi 2011-2019 data without SSC, run without SSC
    'cWO10-rW'% Cal. without Aus and Swi 2011-2019 data without SSC, run with SSC
    'cWO19-rW'% Cal. with Aus and Swi 2011-2019 data without SSC, run with SSC
    'Observed'
    };

legtextsimple = {
    'cW10-rW'
    'cW10-rWO'
    'cWO10-rWO'
    'cWO10-rW'
    'Observed'};

legtextack = {
    'cW10-rW'
    'cW19-rW'
    'cWO10-rWO'
    'cWO19-rWO'
    'Observed'};

path = 'C:\Users\Fredrik.Lagergren\Documents\BioticMod\IPSEuropeSens\Withorwithout\fixwscal\';
regtext = {
    'S_Swe'
    'Switz'
    'Austria'
    'France'};

salvcut = [
    true% Cal. without Aus and Swi 2011-2019 data with SSC, run with SSC
    true% Cal. with Aus and Swi 2011-2019 data with SSC, run with SSC
    false% Cal. without Aus and Swi 2011-2019 data with SSC, run without SSC
    false% Cal. with Aus and Swi 2011-2019 data with SSC, run without SSC
    false% Cal. without Aus and Swi 2011-2019 data without SSC, run without SSC
    false% Cal. with Aus and Swi 2011-2019 data without SSC, run without SSC
    true% Cal. without Aus and Swi 2011-2019 data without SSC, run with SSC
    true% Cal. with Aus and Swi 2011-2019 data without SSC, run with SSC
    ];
    
% Calibration factor for the different settings
k0_news = [
    0.005124114% Cal. without Aus and Swi 2011-2019 data with SSC, run with SSC
    0.002196057% Cal. with Aus and Swi 2011-2019 data with SSC, run with SSC
    0.005124114% Cal. without Aus and Swi 2011-2019 data with SSC, run without SSC
    0.002196057% Cal. with Aus and Swi 2011-2019 data with SSC, run without SSC
    0.000144188% Cal. without Aus and Swi 2011-2019 data without SSC, run without SSC
    9.841E-06% Cal. with Aus and Swi 2011-2019 data without SSC, run without SSC
    0.000144188% Cal. without Aus and Swi 2011-2019 data without SSC, run with SSC
    9.841E-06% Cal. with Aus and Swi 2011-2019 data without SSC, run with SSC
    ];

% Parameter-combination number for the "best" model with the different settings
paras = [
    108382% Cal. without Aus and Swi 2011-2019 data with SSC, run with SSC
    367624% Cal. with Aus and Swi 2011-2019 data with SSC, run with SSC
    108382% Cal. without Aus and Swi 2011-2019 data with SSC, run without SSC
    367624% Cal. with Aus and Swi 2011-2019 data with SSC, run without SSC
    16504% Cal. without Aus and Swi 2011-2019 data without SSC, run without SSC
    11361% Cal. with Aus and Swi 2011-2019 data without SSC, run without SSC
    16504% Cal. without Aus and Swi 2011-2019 data without SSC, run with SSC
    11361% Cal. with Aus and Swi 2011-2019 data without SSC, run with SSC
    ];
    
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
base_bm = 0.0001;

sankillfrac = 0.25; % Fraction of bark beetles killed in sanitary counter-measures
sanmin = 0.01; % Min fraction of killed spruce biomass for sanitary cutting to take place

weights = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\standweightNy2.txt','\t',1,0);

load('Marine_obs_damage_filled_SDI.mat'); % obs_damage2ny(region,year,type), 31 years (1990-2020),
% 26 regions (1-10 S Swe, 11-18 C Swe, 19 Aus, 20-24 Fr, 25-26 Swi),2 types (1=storm 2=Ips)
obs_damage = NaN(26,60,2);
obs_damage(:,30:60,:) = obs_damage2ny;

timedata = NaN(regions,years-10,9); % ,reg,year,setting 1-8 + obs
timedata(:,:,9) = obs_damage(:,:,2);
timeplotdata = NaN(4,years-10,9); % ,large_reg,year,setting 1-8 + obs


% Run the stand-alone model with the different settings
for setting = 1:8
    para = paras(setting);
    k0_new = k0_news(setting);
    if salvcut(setting)
        load('logfileresultsNy2withSSCfix.mat');
        mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\fixedwscalNy2_withSalvSanCut\mwscal_BNE.out','',1,3);
    else
        load('logfileresultsNy2withoutSSCfix.mat');
        mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\fixedwscalNy2_withoutSalvSanCut\mwscal_BNE.out','',1,3);
    end
    gridcells = size(inputdata,1);
    mwscal = reshape(mwscal,[120 gridcells 12]);
    mwscal = permute(mwscal,[3 1 2]);
    wscalm_current = permute(mean(mwscal(5:7,51:120,:)),[3 2 1]);
    wscalm_previous = permute(mean(mwscal(5:7,50:119,:)),[3 2 1]);
    
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
            (inputdata(:,:,:,year,2) > 0) .* salvcut(setting));

        mortality_insect(:,:,:,year) = min(max_insect_mort * inputdata(:,:,:,year,2), k0_new .* patchPI);
    end

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
        
        timedata(region,:,setting) = mod_damage(region,:);
    end
    
end

% Summarize for large regions/countries
timeplotdata(1,:,:) = mean(timedata(1:10,:,:));%Southern Sweden
%timeplotdata(2,:,:) = mean(timedata(11:18,:,:));%Central Sweden
timeplotdata(2,:,:) = mean(timedata(25:26,:,:));%Switzerland
timeplotdata(3,:,:) = timedata(19,:,:);%Austria
timeplotdata(4,:,:) = mean(timedata(20:24,:,:));%France


colors = [
    209	0	255
    0	73	73
    255	182	219
    50	220	90
    150	20	20
    73	0	146
    173	169	74
    76	217	255
    0	0	0
    ];
colors = colors/255;

styles = {'-'
    '-'
    '--'
    '--'
    '-'
    '-'
    '--'
    '--'
    '-'};

markers = 'sdsdsdsdo';

xstart = 0.06;
ystart = 0.09;
xwidth = 0.92;
yheight = 0.88;

xdata = 1961:2020;


% Timeplots
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
        hold(axel1,'all');
        xlim(axel1,[1990 2020]);
        ylabel(axel1,'% damaged');
        hold(axel1,'all');
    
    for setting = 1:9
        ydata = permute(timeplotdata(large_reg,:,setting),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(setting,:),'LineWidth',2,'LineStyle',styles{setting},...
            'Marker',markers(setting),'MarkerEdgeColor','none','MarkerFaceColor',colors(setting,:));
    end
    legend(axel1,legtext,...
        'FontSize',16,'Location','northwest');
    thefile = sprintf('%s',path,'Timecomp',regtext{large_reg},'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end

% Timeplots simplified
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
        hold(axel1,'all');
        xlim(axel1,[1990 2020]);
        ylabel(axel1,'% damaged');
        hold(axel1,'all');
    
    for setting = 1:2:9
        ydata = permute(timeplotdata(large_reg,:,setting),[2 1]) * 100;
        plot(axel1,xdata, ydata,'Color',colors(setting,:),'LineWidth',2,'LineStyle',styles{setting},...
            'Marker',markers(setting),'MarkerEdgeColor','none','MarkerFaceColor',colors(setting,:));
    end
    legend(axel1,legtextsimple,...
        'FontSize',16,'Location','northwest');
    thefile = sprintf('%s',path,'TimecompSimple',regtext{large_reg},'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end

%Ackumulated data
ack_data = NaN(size(timeplotdata));
for large_reg = 1:4
    for setting = 1:9
        firstobs = true;
        for year = 30:60
            if isfinite(timeplotdata(large_reg,year,9)) && isfinite(timeplotdata(large_reg,year,setting))
                if firstobs
                    ack_data(large_reg,year,setting) = timeplotdata(large_reg,year,setting);
                else
                    ack_data(large_reg,year,setting) = ack_data(large_reg,year-1,setting) + timeplotdata(large_reg,year,setting);
                end
                firstobs = false;
            end
        end
    end
end

% Timeplots ackumulated
for large_reg = 1:4
    
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 4]);
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
        hold(axel1,'all');
        xlim(axel1,[1990 2020]);
        ylabel(axel1,'% damaged');
        hold(axel1,'all');
    
    for setting = [1 2 5 6 9]
        ydata = permute(ack_data(large_reg,:,setting),[2 1]) * 100;
        if setting == 1
            LW = 3.5;
            MS = 9;
        else
            LW = 2;
            MS = 6;
        end
        plot(axel1,xdata, ydata,'Color',colors(setting,:),'LineWidth',LW,'LineStyle',styles{setting},...
            'Marker',markers(setting),'MarkerEdgeColor','none','MarkerFaceColor',colors(setting,:),...
            'MarkerSize',MS);
    end
    legend(axel1,legtextack,...
        'FontSize',16,'Location','northwest');
    thefile = sprintf('%s',path,'Ackumulated_bold',regtext{large_reg},'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end

