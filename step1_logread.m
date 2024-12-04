% Reads the patch level output written to the logfile and sort it up in
% Matlab matrice files

clear
close all hidden

%logfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\guesslog20240216.txt'; % DefaultWith
logfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\guesslog20240217.txt'; % DefaultWithout
%logfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\guesslog20240224.txt'; % CalibratedWith
%logfile = 'C:\GuessLibs\crop_ncep_biotics\build\Debug\guesslog20240229.txt'; % CalibratedWith_2plus

% Labels used in the logfiles to distinguish the different types of output
label1 = 'Insect';
label2 = 'Stormpatch';
label3 = 'Phen';
label4 = 'Wscal';

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
years = 2020-1950;

% Insect output matris
vars = 13;
outdata = NaN(gridcells,stands,patches,ifts,years,vars);
% vars: 1=insect_litter, 2=storm_sub, 3=GC_popI, 4=Patch_popI, 5=GC_popInew,
% 6=insect_stem, 7=mort_insect, 8=drougtI, 9=phenI, 10=patchI,
% 11=gridcellI, 12=R, 13=wscal

% Storm output matris
vars2 = 7;
outdata2 = NaN(gridcells,stands,patches,years,vars2);
% var2: Not used year,stand.id, patch.id, patch.managed, patch.thinned,
% 1=sens_ind, 2=dam_frac, 3=dam_wood, 4=cmass_wood, 5=root_mass, 6=height, 7=root_stab_ind

% Insect phenology output matris
vars3 = 6;
outdata3 = NaN(gridcells,years,vars3);
% vars3: 1=last_flightday, 2=swarm2startday, 3=asp, 4=asp30, 5=phenI, 6=storm_rnd

% Matris for compilation of variables used to run the "stand alone" version of the SBB
% damage module
vars4 = 7;
inputdata = NaN(gridcells,stands,patches,years,vars4);
% vars4: 1=insect_litter, 2=insect_stem, 3=wscal, 4=dam_wood, 5=cmass_wood,
% 6=asp, 7=asp30

% Matris for compilation of modelled storm and SBB damage
vars5 = 3;
damfrac_mod = NaN(gridcells,stands,patches,years,vars5);
% vars5: 1=damaged fraction storm, 2=damaged fraction insect, 3=insect_stem


% Reading of LPJ.GUESS wscal output file
mwscal = dlmread('C:\GuessLibs\crop_ncep_biotics\build\out\CalibratedWith_2plus\mwscal_BNE.out','',1,3);% !!with/without change 2
mwscal = reshape(mwscal,[120 gridcells 12]);
mwscal = permute(mwscal,[3 1 2]);
wscalm_current = permute(mwscal(5:7,51:120,:),[3 2 1]);
wscalm_previous = permute(mwscal(5:7,50:119,:),[3 2 1]);

vars6 = 6;
wscaldata = NaN(gridcells,stands,patches,ifts,years,vars6);
% vars6; 1-3 wscal May-July, 4-6 wscal previos year May-July

% Fill wscal data with gridcell averages to make sure there are no missing values
wscaldata(:,:,:,:,:,1:3) = permute(repmat(wscalm_current,[1 1 1 ifts stands patches]),[1 5 6 4 2 3]);
wscaldata(:,:,:,:,:,4:6) = permute(repmat(wscalm_previous,[1 1 1 ifts stands patches]),[1 5 6 4 2 3]);

% Reading the insect variables
fid = fopen(logfile);
theline = 0;
gridcell = 0;
while ~ feof(fid)
    line_ex = fgetl(fid);
    while isempty(line_ex) && ~ feof(fid)
        line_ex = fgetl(fid);
    end
    startstring = textscan(line_ex,'%s',1);
    if strcmp('Commencing',startstring{1})
        gridcell = gridcell + 1;
    end
    if strcmp(label1,startstring{1})
        theline = theline+1;
        %nextstring = textscan(line_ex,'%s %d %d %d %6.4f %4.3f %7.1f %7.1f %7.1f %5.4f %5.4f %4.3f %4.3f %4.3f %4.3f %4.3f');
        indat = textscan(line_ex,'%s %f %f %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f');
        outdata(gridcell,indat{5}+1,indat{6}+1,indat{7}+1,indat{4}-1950,:) = cell2mat(indat(8:20));
    end
    theline = theline+1;
end
fclose(fid);

%Reading the storm variables
fid = fopen(logfile);
theline = 0;
gridcell = 0;
while ~ feof(fid)
    line_ex = fgetl(fid);
    while isempty(line_ex) && ~ feof(fid)
        line_ex = fgetl(fid);
    end
    startstring = textscan(line_ex,'%s',1);
    if strcmp('Commencing',startstring{1})
        gridcell = gridcell + 1;
    end
    if strcmp(label2,startstring{1})
        theline = theline+1;
        indat = textscan(line_ex,'%s %d %d %d %d %d %f %f %f %f %f %f %f');
        outdata2(gridcell,indat{3}+1,indat{4}+1,indat{2}-1950,:) = cell2mat(indat(7:13));
        %NaN(gridcells,stands,patches,years,vars2);
    end
    theline = theline+1;
end
fclose(fid);

%Reading the phenology variables
fid = fopen(logfile);
theline = 0;
gridcell = 0;
while ~ feof(fid)
    line_ex = fgetl(fid);
    while isempty(line_ex) && ~ feof(fid)
        line_ex = fgetl(fid);
    end
    startstring = textscan(line_ex,'%s',1);
    if strcmp('Commencing',startstring{1})
        gridcell = gridcell + 1;
    end
    if strcmp(label3,startstring{1})
        theline = theline+1;
        indat = textscan(line_ex,'%s %d %f %f %f %f %f %f');
        outdata3(gridcell,indat{2}-1950,:) = cell2mat(indat(3:8));
        %NaN(gridcells,stands,patches,years,vars2);
    end
    theline = theline+1;
end
fclose(fid);

%Reading wscal
fid = fopen(logfile);
theline = 0;
gridcell = 0;
while ~ feof(fid)
    line_ex = fgetl(fid);
    while isempty(line_ex) && ~ feof(fid)
        line_ex = fgetl(fid);
    end
    startstring = textscan(line_ex,'%s',1);
    if strcmp('Commencing',startstring{1})
        gridcell = gridcell + 1;
    end
    if strcmp(label4,startstring{1})
        theline = theline+1;
        indat = textscan(line_ex,'%s %f %f %d %d %d %d %f %f %f %f %f %f');
        wscaldata(gridcell,indat{5}+1,indat{6}+1,indat{7}+1,indat{4}-1950,:) = cell2mat(indat(8:13));
        %NaN(gridcells,stands,patches,years,vars2);
    end
    theline = theline+1;
end
fclose(fid);

tempdata = permute(outdata,[1 2 3 5 6 4]);

% Summary of simulated varibles
guessresults = NaN(gridcells,stands,patches,years,7);
guessresults(:,:,:,:,1) = tempdata(:,:,:,:,8);
guessresults(:,:,:,:,2) = tempdata(:,:,:,:,9);
guessresults(:,:,:,:,3) = tempdata(:,:,:,:,3);
guessresults(:,:,:,:,4) = tempdata(:,:,:,:,11);
guessresults(:,:,:,:,5) = tempdata(:,:,:,:,10);
guessresults(:,:,:,:,6) = tempdata(:,:,:,:,4);
guessresults(:,:,:,:,7) = tempdata(:,:,:,:,7) .* tempdata(:,:,:,:,6);

inputdata(:,:,:,:,1) = tempdata(:,:,:,:,1);
inputdata(:,:,:,:,2) = tempdata(:,:,:,:,6);
inputdata(:,:,:,:,3) = tempdata(:,:,:,:,13);
inputdata(:,:,:,:,4:5) = outdata2(:,:,:,:,3:4);
inputdata(:,:,:,:,6:7) = permute(repmat(outdata3(:,:,3:4),[1 1 1 12 5]),[1 4 5 2 3]);
initPI = tempdata(:,:,:,1,4);
initPInext = tempdata(:,:,:,2,4); %Year two data, used if missing for the first year
damfrac_mod(:,:,:,:,1) = outdata2(:,:,:,:,2);
damfrac_mod(:,:,:,:,2) = tempdata(:,:,:,:,7);
damfrac_mod(:,:,:,:,3) = tempdata(:,:,:,:,6);

% Make sure that there are no NaN values
% Insect_litter and insect_stem replaced by 0
tempdata = inputdata(:,:,:,:,1:2);
tempdata(isnan(tempdata)) = 0;
inputdata(:,:,:,:,1:2) = tempdata;
%wscal replaced by 0.95
tempdata = inputdata(:,:,:,:,3);
tempdata(isnan(tempdata)) = 0.95;
inputdata(:,:,:,:,3) = tempdata;
%damaged fraction replaced by 0
damfrac_mod(isnan(damfrac_mod)) = 0;

wscaldata = permute(wscaldata,[1 2 3 5 6 4]); % Take away ift dimension

% Make sure there is a start population
initPI(initPI == 0) = 0.0005; % three decimals in LPJ-GUESS output
initPI(isnan(initPI)) = initPInext(isnan(initPI)); % Replace missing with year two
initPI(isnan(initPI)) = 1; % Should not be needed

%save('logfileresultsNy2callibrated.mat','inputdata','initPI','damfrac_mod','wscaldata');
%save('logfileresultsNy2callibrated2p.mat','inputdata','initPI','damfrac_mod','wscaldata');
%save('logfileresultsNy2withSSCfix.mat','inputdata','initPI','damfrac_mod','wscaldata');
save('logfileresultsNy2withoutSSCfix.mat','inputdata','initPI','damfrac_mod','wscaldata');

%save('guessresultsNy2callibrated.mat','guessresults');
%save('guessresultsNy2callibrated2p.mat','guessresults');
%save('guessresultsNy2withSSCfix.mat','guessresults');
save('guessresultsNy2withoutSSCfix.mat','guessresults');

