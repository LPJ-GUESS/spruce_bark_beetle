% Makes map showing the simulated gridcell with different colour for the
% different parts

clear
close all hidden

addpath('C:\Program Files\MATLAB\R2018b\resources\borders');
% https://climate.copernicus.eu/visualising-and-processing-climate-data-within-matlab

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

gridlist = dlmread('C:\GuessLibs\crop_ncep_biotics\build\Debug\gridlistNy2Clean.txt','\t');

rgbs = [
    219	109	0
    209	0	255
    0	73	73
    255	182	219
    173	169	74
    50	220	90
    150	20	20
    73	0	146
    76	217	255
    0	0	0
    ];
rgbs = rgbs/255;

path = 'C:\Users\Fredrik.Lagergren\Documents\BioticMod\IPSEuropeSens\Withorwithout\';
figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 8 9.8]);

worldmap([43 60],[0 20]);

for region = 1:10
    for gc = cellregion(region,1):cellregion(region,2)
        geoshow(gridlist(gc,2),gridlist(gc,1),'DisplayType','Point','Marker','o','MarkerSize',10,...
            'MarkerFaceColor',rgbs(region,:),'MarkerEdgeColor',rgbs(region,:));
    end
end

for region = 19:26
    for gc = cellregion(region,1):cellregion(region,2)
        geoshow(gridlist(gc,2),gridlist(gc,1),'DisplayType','Point','Marker','o','MarkerSize',10,...
            'MarkerFaceColor',rgbs(region-18,:),'MarkerEdgeColor',rgbs(region-18,:));
    end
end

setm(gca, 'FontSize', 20);
setm(gca,'GLineWidth',1.0);

borders('countries');

thefile = sprintf('%s',path,'map_gridcells');
saveas(figure1,thefile,'emf');
