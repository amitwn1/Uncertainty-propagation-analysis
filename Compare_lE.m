% This script compares betwee nsimulated and observed lE
clear
clc

addpath('C:\PhD\Codes\EC analysis')

% % Define path and file name of the data
[path,EC_filename,biomet_filename] = DefinePaths();

% % Read data (Both raw and filled data)
[EC_table,biomet_table,lE_filled] = ReadData(path,EC_filename,biomet_filename);

% % Get observed LE (aggregated to daily values)
% Convert datetime data to numeric representations
dates_num = convertTo(EC_table.date,'excel');

% Make a vector of unique date values
dates_num_unique = unique(dates_num);

% Get DOYs of observations
temp = unique(floor(EC_table.DOY));

lE_obs_DOY = temp(2:end-1);

% Preallocate the vector containing the daily lE
lE_filled_daily = zeros(length(dates_num_unique)-2,1);

% Aggregate to daily values
n = 1;

for i = 2:length(dates_num_unique)-1   % The first and last days are excluded because they don't contain a full day's data
    
    % Get indices of the current analyzed data
    ind_temp = find(dates_num == dates_num_unique(i));
    
    % Since at each day, the measurement marked with time=00:00 in fact
    % contain data of the previous day, the first and last values in the
    % indices vector are modified
    ind_temp_corrected = [ind_temp(2:end);ind_temp(end)+1];

    % ...
    lE_filled_daily(n) = mean(lE_filled(ind_temp_corrected));   % W/m^2

    n = n + 1;

end

if ~exist('PF_results_summ')
    load('C:\PhD\Codes\PF\particle_simulated_data_240310_16_HAVA1902\final_results_1.mat','PF_results_summ','weights_summ')
end

%% % Get simulated data by DSSAT alone

% lE from DSSAT simulation
lE_DSSAT_weighted_full_season = weights_summ{end}'*PF_results_summ{end}.lE_DSSAT;

% Get the indices of the dates on which lE measurements are available (from
% days of DSSAT simulation)
DSSAT_DOYs = 103:103 + length(lE_DSSAT_weighted_full_season)-1;  % Planting DOY of Gadash 2019 seaons is 103

inds_lE_DSSAT = logical(sum(DSSAT_DOYs == lE_obs_DOY,1));

lE_DSSAT_weighted = lE_DSSAT_weighted_full_season(inds_lE_DSSAT);


%% % Get simulated data by DSSAT-SCOPE (SCOPE calculation)

% % Get weather for SCOPE daily run
% Define weather file name
wth_filename = 'wth_Gadash2019_365_filled.xlsx';

% Get the weather data
[wth_data] = wth_processing_SCOPE_hour(wth_filename);

b.LAI = weights_summ{end}'*PF_results_summ{end}.LAI;  % Take the results from the last DA event
b.SLA = weights_summ{end}'*PF_results_summ{end}.SLAD;
b.LN_conc = weights_summ{end}'*PF_results_summ{end}.leaf_N;

% Preallcations
[lE_hourly_sim,Ac_hourly_sim] = deal(zeros(size(PF_results_summ{1}.LAI,2),24));
[lE_daily_sim,Ac_daily_sim] = deal(zeros(size(PF_results_summ{1}.LAI,2),1));

for i = 1:size(PF_results_summ{1}.LAI,2)
    
    [flux_table] = run_SCOPE_daily_flux_2024(b,i,wth_data,wth_data.z);

    lE_hourly_sim(i,:) = flux_table.lEtot(2:end)';
    Ac_hourly_sim(i,:) = flux_table.Actot(2:end)';

    lE_daily_sim(i) = mean(flux_table.lEtot(2:end),'omitnan');
    Ac_daily_sim(i) = mean(flux_table.Actot(2:end),'omitnan');

end

save('fluxes_SCOPE','lE_hourly_sim','Ac_hourly_sim','lE_daily_sim','Ac_daily_sim')


%% Plot

% Get DOYs of DA
DAPs_of_DA = [13,18,28,33,38,43,58,73,93,98,103,108]; % For Gadash 2019 season

DOYs_of_DA = [13,18,28,33,38,43,58,73,93,98,103,108]+103;  % Planting DOY of Gadash 2019 seaons is 103

planting_DOY = 103;  % For Gadash 2019 season

figure

% Plot observation
plot(lE_obs_DOY-planting_DOY,lE_filled_daily)

hold on

% Plot DSSAT (the unit conversion from mm/day to W/m^2 was already preformed in the PF code)
plot(DSSAT_DOYs - planting_DOY,lE_DSSAT_weighted_full_season)

% Plot SCOPE (coupled model)
plot(0:length(lE_daily_sim)-1,lE_daily_sim);





%% Functions

%% Define paths
function [path,EC_filename,biomet_filename] = DefinePaths()

    % path
    path = 'C:\PhD\Data\Gadash 2019\EC\Processed\';
    % EC data filename (full_output file)
    EC_filename = 'eddypro_2_full_output_2020-08-19T131907_exp';
    % Biomet data filename
    biomet_filename =  'eddypro_2_biomet_2020-08-19T131907_exp.csv';

end

%% Read data
function [EC_table,biomet_table,lE_filled] = ReadData(path,EC_filename,biomet_filename)

% Read EC data file (full output file)
EC_table = importECtable([path,EC_filename]);

% read biomet file
biomet_table = readtable([path,biomet_filename]);

% Read gap filled LE data
filled_Data = readtable('C:\PhD\Codes\EC analysis\lE_processing\Gadash_2019_filled.txt');

lE_filled = filled_Data.LE_f(1:4019);

end