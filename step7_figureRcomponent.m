% Makes figures showing the details of the dependency of water and
% phenology in the SBB module. Presented in the Supplement.

clear
close all hidden

gridfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\gridlistNy2Clean.txt';
% obs_damage(year,region,type), years 1961-2019, 35 regions,2 types 1=storm 2=Ips

gridlist = dlmread(gridfile,'\t');
gridcells = size(gridlist,1);

% Soil moisture data from Copernicus
% https://land.copernicus.eu/global/products/swi
files = {
    'D:\Data\Copernicus\c_gls_SWI-TS_202206300000_C0337_ASCAT_V3.2.1.nc'
    'D:\Data\Copernicus\c_gls_SWI-TS_202206300000_C0355_ASCAT_V3.2.1.nc'
    'D:\Data\Copernicus\c_gls_SWI-TS_202206300000_C0356_ASCAT_V3.2.1.nc'
    'D:\Data\Copernicus\c_gls_SWI-TS_202206300000_C0357_ASCAT_V3.2.1.nc'};

SWIvariable = 'SWI_040';

lons1 = ncread(files{1},'lon');
lons2 = ncread(files{2},'lon');
lons3 = ncread(files{3},'lon');
lons4 = ncread(files{4},'lon');
lats1 = ncread(files{1},'lat');
lats2 = ncread(files{2},'lat');
lats3 = ncread(files{3},'lat');
lats4 = ncread(files{4},'lat');
SWI1 = ncread(files{1},SWIvariable);
SWI2 = ncread(files{2},SWIvariable);
SWI3 = ncread(files{3},SWIvariable);
SWI4 = ncread(files{4},SWIvariable);

% Recalculate lon-lats to corresponding 0.5 degree grid cells
lons1 = round((lons1-0.25)*2,0)/2+0.25;
lons2 = round((lons2-0.25)*2,0)/2+0.25;
lons3 = round((lons3-0.25)*2,0)/2+0.25;
lons4 = round((lons4-0.25)*2,0)/2+0.25;
lats1 = round((lats1-0.25)*2,0)/2+0.25;
lats2 = round((lats2-0.25)*2,0)/2+0.25;
lats3 = round((lats3-0.25)*2,0)/2+0.25;
lats4 = round((lats4-0.25)*2,0)/2+0.25;

days = size(SWI1,1);
SWIgl = NaN(days,gridcells);
inc = NaN(gridcells,1);

liminclude = 0.5; % Fraction of daily OK values needed to accept gridcell

% Loop over LPJ-GUESS gridcells to calculate average SWI in those cells
for gridcell = 1:gridcells
    lon = gridlist(gridcell,1);
    lat = gridlist(gridcell,2);
    tempsum = zeros(days,1);
    counter = 0;
    for rownr = 1:size(lons1)
        if lon == lons1(rownr) && lat == lats1(rownr) && sum(isfinite(SWI1(:,rownr)))>days*liminclude
            tempsum = tempsum + SWI1(:,rownr);
            counter = 1 + counter;
        end
    end
    for rownr = 1:size(lons2)
        if lon == lons2(rownr) && lat == lats2(rownr) && sum(isfinite(SWI2(:,rownr)))>days*liminclude
            tempsum = tempsum + SWI2(:,rownr);
            counter = 1 + counter;
        end
    end
    for rownr = 1:size(lons3)
        if lon == lons3(rownr) && lat == lats3(rownr) && sum(isfinite(SWI3(:,rownr)))>days*liminclude
            tempsum = tempsum + SWI3(:,rownr);
            counter = 1 + counter;
        end
    end
    for rownr = 1:size(lons4)
        if lon == lons4(rownr) && lat == lats4(rownr) && sum(isfinite(SWI4(:,rownr)))>days*liminclude
            tempsum = tempsum + SWI4(:,rownr);
            counter = 1 + counter;
       end
    end
    SWIgl(:,gridcell) = tempsum ./ counter;
    inc(gridcell) = counter;
end

weights = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\standweightNy2.txt','\t',1,0);


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

load('logfileresultsNy2withSSC.mat');
% inputdata(gridcell,stand,patch,year,var)
% vars: 1=insect_litter, 2=insect_stem, 3=wscal, 4=dam_wood, 5=cmass_wood,
% 6=asp, 7=asp30
% initPI(gridcell,stand,patch)

mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\r2020Ny2_withSalvSanCut\mwscal_BNE.out','',1,3);
mwscal = reshape(mwscal,[120 231 12]);
mwscal = permute(mwscal,[3 1 2]);
wscalm_current = permute(mean(mwscal(5:7,51:120,:)),[3 2 1]);
wscalm_previous = permute(mean(mwscal(5:7,50:119,:)),[3 2 1]);
wscal = ((wscalm_current + wscalm_previous * 3) / 4)';

mwscal = reshape(mwscal,[120*12 231]);
mwscal_reg = NaN(120*12,regions);

asp = inputdata(:,:,:,:,6);
asp = permute(mean(mean(permute(asp,[2 3 4 1]))),[3 4 1 2]);
asp30 = inputdata(:,:,:,:,7);
asp30 = permute(mean(mean(permute(asp30,[2 3 4 1]))),[3 4 1 2]);
phenR = asp * 3.33 ./ (asp30 + 75);
wscalR_5 = min(2,(1-wscal) * 5);
wscalR_10 = min(2,(1-wscal) * 10);
wscalR_27 = min(2,(1-wscal) * 27);

load('Marine_obs_damage_filled_SDI.mat'); % obs_damage2ny(region,year,type), 31 years (1990-2020),
% 26 regions (1-10 S Swe, 11-18 C Swe, 19 Aus, 20-24 Fr, 25-26 Swi),2 types (1=storm 2=Ips)
obs_damage = NaN(26,60,2);
obs_damage(:,30:60,:) = obs_damage2ny;

years = size(asp,1);

asp_reg = NaN(years,regions);
phenR_reg = NaN(years,regions);
wscal_reg = NaN(years,regions);
wscalR5_reg = NaN(years,regions);
wscalR10_reg = NaN(years,regions);
wscalR27_reg = NaN(years,regions);
SWI_reg = NaN(days,regions);

% Calculate averages for regions weighted by forest fraction in gridcell
for region = 1:regions
    startcell = cellregion(region,1);
    endcell = cellregion(region,2);
    cells = endcell - startcell + 1;
    gc_weights = weights(startcell:endcell,13); % Column 13 is the total fraction of forest
    gc_weights = gc_weights ./ sum(gc_weights);
    asp_reg(:,region) = sum(asp(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    phenR_reg(:,region) = sum(phenR(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    wscal_reg(:,region) = sum(wscal(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    wscalR5_reg(:,region) = sum(wscalR_5(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    wscalR10_reg(:,region) = sum(wscalR_10(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    wscalR27_reg(:,region) = sum(wscalR_27(:,startcell:endcell) .* repmat(gc_weights',[years 1]),2);
    SWI_reg(:,region) = sum(SWIgl(:,startcell:endcell) .* repmat(gc_weights',[days 1]),2);
    mwscal_reg(:,region) = sum(mwscal(:,startcell:endcell) .* repmat(gc_weights',[120*12 1]),2);
end

damfrac = permute(obs_damage(:,:,2),[2 1]);
invR = log(damfrac(2:60,:) ./ damfrac(1:59,:));

xstart = 0.09;
xwidth = 0.84;
yheight = 0.29;

xyears = (1980:2020)';
xmonths = ((1:41*12)/12 + 1980 - 0.5)';
xdays = ((1:365*41+11)/365.2683 + 1980 - 0.5)';

path = 'C:\Users\Fredrik.Lagergren\Documents\BioticMod\IPSEuropeSens\Withorwithout\';
colors = [0.9 0.5 0.1
    0.1 0.25 1
    0.1 0.3 0.1
    0.8 0.4 0.8];
markers = 'osd^';


for region = 1:regions
    figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 12 16]);
    
    ystart = 0.05;
    
    % Bottom graphs left y-axis: Inventoried SBB damaged fraction
    axel3 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
    xlim(axel3,[1999 2021]);
    hold(axel3,'all');
    ydata = damfrac(20:60,region) * 100;
    yyaxis(axel3,'left');
    plot(axel3,xyears, ydata,'Color',colors(4,:),'LineWidth',2,'Marker','o',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(4,:),'MarkerSize',10);
    ylabel(axel3,'Damaged fraction (%)');
    hold(axel3,'all');
    
    % Bottom graphs right y-axis: SBB increase rate calculated from damaged fraction
    yyaxis(axel3,'right');
    ydata = invR(19:59,region);
    plot(axel3,xyears, ydata,'Color',colors(1,:),'LineWidth',2,'Marker','^',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerSize',10);
    ylim(axel3,[-2.5 2.5]);
    ylabel(axel3,'Increase rate (R)');
    xlabel(axel3,'Year');
    hold(axel3,'all');
    ax = gca;
    ax.YAxis(1).Color = colors(4,:);
    ax.YAxis(2).Color = colors(1,:);
    hold(axel3,'all');
    
    ystart = 0.37;
    
    % Middle graphs: Simulated wscal (Rwater, with three settings of k5) and
    % phenology (Rphen) depending components of SBB population index grow rate
    axel1 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
    xlim(axel1,[1999 2021]);
    hold(axel1,'all');
    ydata = phenR_reg(30:70,region);
    plot(axel1,xyears, ydata,'Color',colors(1,:),'LineWidth',2,'Marker','o',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerSize',10);
    hold(axel1,'all');
    ydata = wscalR5_reg(30:70,region);
    plot(axel1,xyears, ydata,'Color',colors(2,:),'LineWidth',2,'Marker','s',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerSize',10);
    hold(axel1,'all');
    ydata = wscalR10_reg(30:70,region);
    plot(axel1,xyears, ydata,'Color',colors(2,:),'LineWidth',2,'Marker','d',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerSize',10);
    ydata = wscalR27_reg(30:70,region);
    plot(axel1,xyears, ydata,'Color',colors(2,:),'LineWidth',2,'Marker','^',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerSize',10);
    ylim(axel1,[0 2]);
    ylabel(axel1,'Increase rate (R) component');
    hold(axel1,'all');
    
    ystart = 0.69;
    
    % Top graphs left y-axis: Daily SWI and wscal monthly and averaged over
    % summer month for previous and current year
    axel2 = axes('Position',[xstart ystart xwidth yheight],'Box','on','XGrid','on','YGrid','on',...
        'FontSize',16,'Parent',figure1);
    xlim(axel2,[1999 2021]);
    hold(axel2,'all');
    ydata = NaN(size(xdays));
    yyaxis(axel2,'left');
    ydata(9863:14976) = SWI_reg(1:5114,region) *0.01;
    plot(axel2,xdays, ydata,'Color',colors(3,:),'LineWidth',2,'LineStyle','-','Marker','none',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(3,:));
    hold(axel2,'all');
    ydata = mwscal_reg(79*12+1:1440,region);
    plot(axel2,xmonths, ydata,'Color',colors(2,:),'LineWidth',2,'LineStyle','-','Marker','o',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerSize',5);
    hold(axel2,'all');
    ydata = wscal_reg(30:70,region);
    plot(axel2,xyears, ydata,'Color',colors(2,:),'LineWidth',2,'LineStyle','-','Marker','^',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(2,:),'MarkerSize',10);
    hold(axel2,'all');
    ylim(axel2,[0 1]);
    ylabel(axel2,'Wscal, SWI');
    hold(axel2,'all');
    
    % Top graphs right y-axis: Length of Autumn Swarm Period
    yyaxis(axel2,'right');
    ydata = asp_reg(30:70,region);
    plot(axel2,xyears, ydata,'Color',colors(1,:),'LineWidth',2,'Marker','d',...
        'MarkerEdgeColor','none','MarkerFaceColor',colors(1,:),'MarkerSize',10);
    ylim(axel2,[0 90]);
    ylabel(axel2,'Autumn swarm period (days)');
    hold(axel2,'all');
    ax = gca;
    ax.YAxis(1).Color = colors(2,:);
    ax.YAxis(2).Color = colors(1,:);
    hold(axel3,'all');
    
    thefile = sprintf('%s',path,'Components',num2str(region),'.emf');
    saveas(figure1,thefile,'emf');
    delete(figure1);
    
end


