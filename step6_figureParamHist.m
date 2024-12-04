% Making histogram figure showing the number of the top 50 models, with the best combined
% statistics of bias, RMSE and R2, diveded in order number of the tested
% parameters. Presented in the Supplement.

clear
close all hidden


%The number of the top 50 models, with the best combined statistics of bias, RMSE and R2, diveded in order number (row) of the tested parameters
%Data from step3 exported from Excel

%With Salvage and Sanetary Cutting
parcountW = [
%   All four, cal. data -2010   S Swe, cal. data -2010      Switz, cal. data -2010      Austria, cal. data -2010    NE France, cal. data -2010  All four, cal. data -2019   S Swe, cal. data -2019      Switz, cal. data -2019      Austria, cal. data -2019    NE France, cal. data -2019   
%   k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5
    44	6	0	11	2	0	28	7	50	0	13	5	0	0	5	47	3	1	8	1	3	9	45	35	0	0	0	4	35	0	44	0	0	0	0	3	50	0	47	6	0	0	7	50	0	2	10	0	0	7	50	0	0	6	0	0	23	0	40	0	0	0	0	39	0	47	0	0	1	0
    0	0	0	12	1	0	15	3	0	0	13	7	0	0	5	0	1	2	6	0	1	23	0	13	0	0	3	1	12	0	6	0	0	3	0	6	0	0	3	6	0	0	6	0	0	11	12	0	0	6	0	0	0	8	0	0	16	50	10	0	3	0	0	8	0	3	0	0	5	0
    0	0	0	12	1	0	1	8	0	10	13	7	0	0	7	0	2	0	8	0	0	3	0	2	0	0	3	0	3	0	0	0	0	4	0	7	0	0	0	8	0	0	8	0	8	16	5	0	0	10	0	0	6	6	0	0	7	0	0	0	3	1	0	3	0	0	0	0	9	0
    0	4	15	12	0	0	0	5	0	18	11	13	0	0	12	0	17	19	7	0	0	7	0	0	45	8	3	0	0	39	0	0	2	10	0	10	0	0	0	8	22	0	9	0	17	19	11	0	0	7	0	4	40	9	0	0	2	0	0	50	5	3	0	0	7	0	0	2	11	0
    1	6	13	1	1	2	0	4	0	12	0	13	0	0	7	3	19	26	6	22	46	5	0	0	0	10	8	0	0	11	0	0	6	12	0	11	0	0	0	7	28	22	5	0	16	2	4	0	0	8	0	18	4	9	0	50	2	0	0	0	9	7	0	0	43	0	3	5	13	0
    1	13	5	1	7	18	6	12	0	10	0	4	3	0	7	0	8	1	7	13	0	3	0	0	2	13	16	0	0	0	0	3	14	11	50	7	0	17	0	8	0	28	7	0	9	0	5	0	0	7	0	23	0	7	17	0	0	0	0	0	10	16	0	0	0	0	26	10	8	50
    4	21	17	1	38	30	0	11	0	0	0	1	47	50	7	0	0	1	8	14	0	0	5	0	3	19	17	45	0	0	0	47	28	10	0	6	0	33	0	7	0	0	8	0	0	0	3	50	50	5	0	5	0	5	33	0	0	0	0	0	20	23	50	0	0	0	21	33	3	0
    ];

%Without Salvage and Sanetary Cutting
parcountWO = [
%   All four, cal. data -2010   S Swe, cal. data -2010      Switz, cal. data -2010      Austria, cal. data -2010    NE France, cal. data -2010  All four, cal. data -2019   S Swe, cal. data -2019      Switz, cal. data -2019      Austria, cal. data -2019    NE France, cal. data -2019   
%   k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5  k1  k2  kgc pyw k3  k4  k5
    8	35	0	0	3	0	0	20	50	1	1	7	0	1	15	37	0	17	1	0	0	0	0	0	0	11	0	0	13	0	0	0	5	0	0	16	50	0	0	13	0	0	20	50	6	3	6	0	6	12	42	0	38	4	11	0	6	0	28	2	0	0	0	13	0	0	0	2	0	0
    10	15	0	0	5	0	0	16	0	0	0	10	0	0	13	4	0	7	2	6	0	0	0	0	0	16	0	0	7	0	0	0	7	0	0	17	0	0	0	8	0	0	14	0	0	3	9	0	0	10	5	0	12	5	14	0	6	46	14	0	0	0	0	14	0	0	0	6	0	0
    8	0	0	0	8	0	0	6	0	0	6	6	0	0	11	4	0	1	3	6	0	0	8	0	0	13	0	0	17	0	0	0	9	0	0	6	0	1	0	5	0	0	7	0	0	2	7	0	0	9	3	0	0	4	14	0	12	4	8	0	1	5	0	15	0	0	0	8	0	0
    10	0	0	0	8	0	0	8	0	5	12	3	0	0	7	1	0	11	8	9	0	0	14	0	25	10	0	0	12	0	0	0	9	0	0	10	0	5	0	7	0	0	7	0	7	11	1	0	0	6	0	2	0	7	7	0	8	0	0	3	5	7	0	7	0	0	0	10	0	0
    8	0	0	1	10	31	29	0	0	30	16	6	0	0	3	0	9	5	7	20	0	14	14	0	21	0	50	0	1	1	3	6	11	21	0	0	0	25	6	7	0	0	0	0	27	19	8	0	0	5	0	13	0	9	4	0	8	0	0	1	10	10	0	1	0	22	4	11	50	0
    4	0	0	13	7	19	21	0	0	14	8	9	8	0	1	2	26	7	15	8	21	18	8	5	4	0	0	50	0	15	12	14	7	29	50	0	0	17	21	4	21	0	1	0	10	7	13	7	0	4	0	19	0	9	0	0	8	0	0	2	14	15	4	0	5	25	17	8	0	50
    2	0	50	36	9	0	0	0	0	0	7	9	42	49	0	2	15	2	14	1	29	18	6	45	0	0	0	0	0	34	35	30	2	0	0	1	0	2	23	6	29	50	1	0	0	5	6	43	44	4	0	16	0	12	0	50	2	0	0	42	20	13	46	0	45	3	29	5	0	0
    ];

parcount = NaN(7,7,5,4); % (parvalnr,parameter,region,cal_setting)
parcount(:,:,:,1:2) = reshape(parcountW,[7 7 5 2]);
parcount(:,:,:,3:4) = reshape(parcountWO,[7 7 5 2]);

rgbs = [
    1	1	1
    1	0.71	0.86
    0.43	0.71	1
    0.71	0.43	1
    0.57	0.29	0
    0	0.29	0.29
    0	0	0
    ];

xtotal = 0.96;
xgap = 0.015;
xwidth = (xtotal - xgap * 3) / 4;
ytotal = 0.985;
ystart = 0.025;
yheight = ytotal/7;
ygap = yheight/5 * 1/5;
ysubheight = (yheight - ygap * 9) / 5;


figure1 = figure('Color',[1 1 1],'PaperPosition',[0.6345 0.6345 14 18]);
subfig = 0;
for parameter = 1:7
    
    for region = 1:5
        ysubstart = ystart + (5 - region) * (ysubheight + ygap);
        xstart = 0.02;
        
        for cal_setting = 1:4
            subfig = subfig + 1;
            axel.subfig = axes('Position',[xstart ysubstart xwidth ysubheight],...
                'Box','on','XGrid','on','YGrid','on','FontSize',16,'Parent',figure1);
            bar(axel.subfig,parcount(:,parameter,region,cal_setting),'FaceColor',rgbs(region+1,:));
            xlim(axel.subfig,[0.4 7.6]);
            ylim(axel.subfig,[0 50]);
            hold(axel.subfig,'all');
            
            axel.subfig.YTickLabel = ' ';
            if region < 5
                axel.subfig.XTickLabel = ' ';
            end
            axel.subfig.FontSize = 20;
            hold(axel.subfig,'all');
            
            xstart = xstart + xwidth + xgap;
        end
        
    end
    ystart = ystart + yheight;
    
end
        
saveas(figure1,'C:\Users\Fredrik.Lagergren\Documents\BioticMod\manus\Testing\parcountfig.emf','emf');
