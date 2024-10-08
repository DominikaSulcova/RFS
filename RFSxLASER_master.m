%% RFS vs. laser - master script
% author:   Dominika Sulcova
%           MSH - Medical School Hamburg
% created:  2024   
% 
% This is a master script for a study comparing psychophisics and
% neurphysiological response to painful stimulation using a gold standard
% YAP-laser stimulation and a novel radiofrequency stimulator (RFS)
% develpped by the engineering faculty of the University of Oro Verde (AR).
% 
% Collected data: 
%   1) evoked potentials 
%           --> 2 stimulators: laser (LEP), RFS (RFEP)
%           --> 2 intensities: threshold (low), 1.5 * threshold (high)
%           --> 2 areas: left hand (left), right hand (right)
%   2) psychophysical evaluation of painful stimulation
%           --> threshold values (determined by pain-related descriptors)
%           --> intensity ratings acquired at single trial level 
%   3) resting state EEG (RS) 
%           --> 1.5 mins eyes open + 1.5 mins eyes closed
%           --> 2 repetitions: before / after stimulation
% 
% The first part of the script runs during the experimental session and
% allows to:    1) encode subject data
%               2) perform threshold tracking 
%               3) calculate the pain threshold based on psychophysical
%               descriptors --> imported from PsychoPy as .csv
% ==> output:   structure RFSxLASER_info
%               within the output file RFSxLASER_output.mat
% 
% The second part of the script performs EEG pre-processing of single
% subject data with follwing steps:
%   1)
% ==> output:   info, parameters -->    RFSxLASER_info 
%               measured variables -->  RFSxLASER_measures
%               processed data -->  	RFSxLASER_data
%               within the output file RFSxLASER_output.mat
% 
% The third part of the script provides group average visualization and
% exports the data in appropriate formates for the stat analysis in R.

% ===================== PART 1: experimental session =====================
%% params  
clc, clear all

% directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');    % toolboxes
folder.raw = uigetdir(pwd, 'Coose the input folder');           % raw data --> at MSH, this should be the study folder at the V drive
folder.processed = uigetdir(pwd, 'Choose the data folder');     % processed data --> wherever you want to store the voluminous EEG data
folder.output = uigetdir(pwd, 'Choose the OneDrive folder');    % output folder --> figures, loutput file, exports 
cd(folder.output)

% output
study = 'RFSxLASER';
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
figure_counter = 1;

% sound
load handel.mat
haleluja = y; clear y Fs

%% 1) fill in subject & session info 
% ----- section input -----
params.laser.pulse = 3;
params.laser.diameter = 5;
params.RFS.pulse = 10;
params.RFS.probe = 5;
% -------------------------

% load the info structure
if exist(output_file) == 2
    output_vars = who('-file', fileName);
    if ismember('RFSxLASER_info', output_vars)
        load(output_file, 'RFSxLASER_info')
    else
        RFSxLASER_info = struct;
        save(output_file, 'RFSxLASER_info','-append')
    end
else
    RFSxLASER_info = struct;
    save(output_file, 'RFSxLASER_info')
end

% get subject & session info
prompt = {'date:', 'subject:', 'age:', 'male:', 'handedness:', 'height (cm):', 'weight (kg):', 'right arm (cm):', 'left arm (cm):'};
dlgtitle = 'Subject session information';
dims = [1 35];
definput = {date, 'S000', '18', '1', '1', '175', '70', '60', '60'};
session_info = inputdlg(prompt,dlgtitle,dims,definput);
clear prompt dlgtitle dims definput

% identify subject's index
subject_idx = str2num(session_info{2}(end-1:end));

% fill in the metadata structure
RFSxLASER_info(subject_idx).date = session_info{1};
RFSxLASER_info(subject_idx).ID = session_info{2};
RFSxLASER_info.age = str2num(session_info{3});
RFSxLASER_info.male = str2num(session_info{4});
RFSxLASER_info.handedness = str2num(session_info{5});
RFSxLASER_info(subject_idx).body.height = str2num(session_info{6});
RFSxLASER_info(subject_idx).body.weight = str2num(session_info{7});
RFSxLASER_info(subject_idx).body.arm_right = str2num(session_info{8});
RFSxLASER_info(subject_idx).body.arm_left = str2num(session_info{9});
RFSxLASER_info(subject_idx).stimulation.laser.pulse = params.laser.pulse;
RFSxLASER_info(subject_idx).stimulation.laser.diameter = params.laser.diameter;
RFSxLASER_info(subject_idx).stimulation.RFS.pulse = params.RFS.pulse;
RFSxLASER_info(subject_idx).stimulation.RFS.probe = params.RFS.probe;
clear session_info

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');
clear params

%% 2) run threshold tracking 
% ----- section input -----
params.stimulus = {'laser', 'RFS'};
% -------------------------

% run for both stimuli
for s = 1:length(params.stimulus)
    % press SPACE to start timing the session 
    fprintf('Threshold tracking: %s\n', params.stimulus{s});
    ok = 0; 
    while ok == 0    
        key = input('Press ENTER to start', 's');
        if isempty(key)  
            ok = 1;  
        end
    end    
    fprintf('Let''s go!');
    clear ok key  

    % initiate variables
    intensity = struct;
    trial_counter = 1;
    continue_trials = true;
    
    % threshold tracking - loop through trials
    while continue_trials
        % encode current intensity
        prompt = {'intensity:'};
        dlgtitle = sprintf('trial n. %d', trial_counter);
        dims = [1 35];
        definput = {'00'};
        trial.intensity = (inputdlg(prompt,dlgtitle,dims,definput));
        clear prompt dlgtitle dims definput
        
        % determine time
        trial.time = sprintf('%s', datetime);
        
        % append to the output vector
        intensity(trial_counter).trial = trial_counter;
        intensity(trial_counter).time = trial.time(end-7:end);
        intensity(trial_counter).intensity = str2num(trial.intensity{:});
      
        % ask to continue
        answer = questdlg('Continue threshold tracking?', 'next step',...
            'Continue', 'Finish', 'Continue');
        if strcmp(answer, 'Continue')
            % update trial counter
            trial_counter = trial_counter + 1;
        else
            % leave the loop
            continue_trials = false;
        end
    end

    % append to info structure
    if strcmp(params.stimulus{s}, 'laser')
        RFSxLASER_info(subject_idx).stimulation.laser.threshold.unit = 'J';
        RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures = intensity;
    elseif strcmp(params.stimulus{s}, 'RFS')
        RFSxLASER_info(subject_idx).stimulation.RFS.threshold.unit = '%MSO';
        RFSxLASER_info(subject_idx).stimulation.RFS.threshold.measures = intensity;
    end
end

% save to the output file
save(output_file, 'RFSxLASER_info', '-append');
clear params s intensity trial_counter continue_trials trial

%% 3) determine thresholds 
% ----- section input -----
params.stimulus = {'laser', 'RFS'};
params.labels = {'nothing' 'warm' 'pricking' 'burning' 'electric'};
params.values = [0 1 2 2 3];
% -------------------------

% append the descriptors
RFSxLASER_info(subject_idx).stimulation.descriptors.labels = params.labels;
RFSxLASER_info(subject_idx).stimulation.descriptors.values = params.values;

% run for both stimuli
for a = 1:length(params.stimulus)
    % provide update
    fprintf('%s pain threshold:\n', params.stimulus{a})

    % import PsychoPy .csv file
    fprintf('importing PsychoPy data...\n')
    file2import = dir(sprintf('%s\\%s\\threshold\\*%s*.csv', folder.raw, RFSxLASER_info(subject_idx).ID, params.stimulus{a}));
    if size(file2import, 1) > 1
        fprintf('ERROR: multiple %s thresholding files attributed to this subject!\n', params.stimulus{a})
        exit()
    end
    option = detectImportOptions(sprintf('%s\\%s', file2import.folder, file2import.name));
    option.SelectedVariableNames = {'trial', 'fixation_stopped', 'descriptor', 'pain'};
    ratings_table = readtable(sprintf('%s\\%s', file2import.folder, file2import.name), option);
    clear file2import option

    % adjust ratings format
    fprintf('extracting trials...\n')
    for b = 1:height(ratings_table)
        % if it is a valid trial
        if ~isnan(ratings_table{b, 1})
            % encode trial
            ratings(b).trial = ratings_table{b, 1};
            
            % encode timestamp
            ratings(b).time_psychopy = ratings_table{b, 2};
            
            % encode descriptors
            desc = ratings_table{b, 3};
            desc = regexprep(desc{:}, '[[]'',]' , '');
            ratings(b).descriptors = split(desc, ' ');
            
            % encode pain intensity
            int = ratings_table{b, 4};
            int = regexprep(int{:}, '[[]'',]' , '');
            ratings(b).pain = split(int, ' ');
        end
    end

    % check if the trial number matches
    if length(RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures) ~= length(ratings)
        fprintf('ERROR: PsychoPy trials do not match MATLAB trials!\n')
        exit()
    end

    % clasify trials and append to output structure
    fprintf('%d valid trials found\nclassifying ...\n', length(ratings)) 
    for c = 1:length(ratings)
        % if trial number matches
        if RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures(c).trial == ratings(c).trial
            % append 
            RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures(c).time_psychopy = ratings(c).time_psychopy;
            RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures(c).descriptors = ratings(c).descriptors;
            RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures(c).pain = ratings(c).pain;
            
            % determine level of evoked sensation based on descriptors
            for d = 1:length(ratings(c).descriptors)
                index = find(strcmp(ratings(c).descriptors{d}, params.labels));
                if ~isempty(index)
                    desc_int(d) = params.values(index);
                else
                    disp('Error: Trial %d contains an unknown descriptor.');
                    exit()
                end                
            end
            RFSxLASER_info(subject_idx).stimulation.laser.threshold.measures(c).intensity = max(desc_int);
            
            % prepare visualization input
            visual.sensation(c) = max(desc_int);
            if visual.intensity(c) == 0
                visual.color(c, :) = [0.6510    0.6510    0.6510];
            elseif visual.intensity(c) == 1
                visual.color(c, :) = [0    0.4471    0.7412];
            elseif visual.intensity(c) == 2
                visual.color(c, :) = [0.8000    0.0157    0.0157];
            elseif visual.intensity(c) == 3
                visual.color(c, :) = [0.8000    0.0157    0.0157];
            end  
            visual.intensity(c) = output_info(subject_idx).psychophysics.threshold_trials(c).power;    
            
            % save to the output file
            save(output_file, 'RFSxLASER_info', '-append');
        end

        
    end

end

clear a b c d file2import ratings_table desc int index desc_int






% identify intenisity transitions
threshold.perception = [];
threshold.pain = [];
for d = 2:length(ratings)
    if intensity(d) == 1 & intensity(d-1) == 0
        threshold.perception(end+1) = mean(power([d-1, d]));
    elseif intensity(d) == 0 & intensity(d-1) == 1
        threshold.perception(end+1) = mean(power([d-1, d]));
    elseif intensity(d) == 2 & intensity(d-1) == 1
        threshold.pain(end+1) = mean(power([d-1, d]));
    elseif intensity(d) == 1 & intensity(d-1) == 2
        threshold.pain(end+1) = mean(power([d-1, d]));
    end
end
output_info(subject_idx).psychophysics.threshold_perception = round(mean(threshold.perception));
output_info(subject_idx).psychophysics.threshold_pain = round(mean(threshold.pain));

% clasify trials based on descriptors
for c = 1:length(ratings)
    % if trial number matches
    if output_info(subject_idx).psychophysics.threshold_trials(c).trial == ratings(c).trial
        % determine value of descriptors 
        int_descriptors = [];
        for d = 1:length(descriptors.labels)             
            if ismember(descriptors.labels(d), output_info(subject_idx).psychophysics.threshold_trials(c).descriptors) 
                int_descriptors(end+1) = descriptors.values(d);
            end
        end
        
        % determine overall intensity 
        output_info(subject_idx).psychophysics.threshold_trials(c).intensity_descriptors = max(int_descriptors);
    end
end

% identify intenisity transitions
intensity_descriptors = extractfield(output_info(subject_idx).psychophysics.threshold_trials, 'intensity_descriptors');
threshold.perception = [];
threshold.pain = [];
for d = 2:length(intensity_descriptors)
    if intensity_descriptors(d) == 1 && intensity_descriptors(d-1) == 0
        threshold.perception(end+1) = mean(power([d-1, d]));
    elseif intensity_descriptors(d) == 0 && intensity_descriptors(d-1) == 1
        threshold.perception(end+1) = mean(power([d-1, d]));
    elseif intensity_descriptors(d) == 2 && intensity_descriptors(d-1) == 1
        threshold.pain(end+1) = mean(power([d-1, d]));
    elseif intensity_descriptors(d) == 1 && intensity_descriptors(d-1) == 2
        threshold.pain(end+1) = mean(power([d-1, d]));
    end
end
output_info(subject_idx).psychophysics.threshold_perception_descriptors = round(mean(threshold.perception));
output_info(subject_idx).psychophysics.threshold_pain_descriptors = round(mean(threshold.pain));

% launch the figure
x = 1:length(ratings);
fig = figure(figure_counter);
set(fig, 'units','normalized','outerposition',[0 0 1 1])

% plot the line
plot(x, power, ':', 'color', [0.6510    0.6510    0.6510], 'linewidth', 3)
hold on

% plot the trials
for e = 1:length(ratings)
    plot(x(e), power(e), 'o', 'MarkerSize', 12, 'MarkerEdgeColor','none','MarkerFaceColor', marker_color(e, :))
    hold on
end

% plot thresholds
line(x(1)-5:x(end)+1, repelem(output_info(subject_idx).psychophysics.threshold_perception, 1, length(x)+6), 'color', 'black', 'linewidth', 1.5)
line(x(1)-5:x(end)+1, repelem(output_info(subject_idx).psychophysics.threshold_pain, 1, length(x)+6), 'color', 'black', 'linewidth', 1.5)
line(x(1)-5:x(end)+1, repelem(output_info(subject_idx).psychophysics.threshold_perception_descriptors, 1, length(x)+6), 'color', [0.6353    0.0784    0.1843], 'linewidth', 1.5)
line(x(1)-5:x(end)+1, repelem(output_info(subject_idx).psychophysics.threshold_pain_descriptors, 1, length(x)+6), 'color', [0.6353    0.0784    0.1843], 'linewidth', 1.5)

% plot annotations
text(x(1)-5 + 0.25, output_info(subject_idx).psychophysics.threshold_perception + 0.25, 'perception threshold', 'fontsize', 14);
text(x(1)-5 + 0.25, output_info(subject_idx).psychophysics.threshold_pain + 0.25, 'pain threshold', 'fontsize', 14);
text(x(1)-5 + 0.25, output_info(subject_idx).psychophysics.threshold_perception_descriptors + 0.25, 'perception threshold - descriptors', 'fontsize', 14, 'color', [0.6353    0.0784    0.1843]);
text(x(1)-5 + 0.25, output_info(subject_idx).psychophysics.threshold_pain_descriptors + 0.25, 'pain threshold - descriptors', 'fontsize', 14, 'color', [0.6353    0.0784    0.1843]);

% adjust visuals
xlim([x(1)-5, x(end)+1]); 
ylim([0, 1.1*max(power)])
xlabel('trials'); ylabel('power (%MSO)')
title(sprintf('threshold tracking - %s', output_info(subject_idx).ID))
set(gca, 'fontsize', 14)
    
% save figure
saveas(fig, sprintf('%s\\figures\\threshold_%s.png', folder.output, output_info(subject_idx).ID))

% update figure counter
figure_counter = figure_counter + 1;
clear descriptors a b c d e int desc ratings intensity power threshold marker_color int_descriptors x fig  

%% block psychophysics


%% ===================== PART 2: single subject processing ================
% params - dataset
clc, clear all
study = 'RFS_OV2024';
figure_counter = 1;

% params - directories
folder.toolbox = uigetdir(pwd, 'Choose the toolbox folder');        % letswave masterfiles
folder.input = uigetdir(pwd, 'Coose the input folder');             % raw data
folder.data = uigetdir(pwd, 'Choose the data folder');              % processed data
folder.output = uigetdir(pwd, 'Choose the output folder');        % output folder --> figures, logfiles, output .mat file
cd(folder.data)

% params - output
output_file = sprintf('%s\\%s_output.mat', folder.output, study);
load(output_file, 'output_info');
% % output_info = struct;
% % save(output_file, 'output_info')

%% import data for letswave
% add letswave 6 to the top of search path
addpath(genpath([folder.toolbox '\letswave 6']));
fprintf('Loading:\n')

% identify the appropriate files
cd(folder.input)
file2import = dir(sprintf('*%s*.txt', output_info(subject_idx).ID));
    
% cycle through datasets
for a = 1:length(file2import)
    % create name for letswave
    data_name = regexprep(file2import(a).name, '-' , '');
    data_name = regexprep(data_name, '_' , ' ');
    data_name(end-19:end) = [];

    % provide update
    fprintf('%s ...\n', data_name)

    % encode the filename to metadata
    output_info(subject_idx).dataset(a).file = file2import(a).name;
    output_info(subject_idx).dataset(a).name = data_name;

    % load the txt file
    data_mat = readmatrix(file2import(a).name);
    data_mat = data_mat(:, 2:9)';
    data_mat_lw = [];
    for b = 1:size(data_mat, 1) - 1
        data_mat_lw(1, b, 1, 1, 1, :) = data_mat(b, :)*(-1);
    end
    data_mat_lw(1, size(data_mat, 1), 1, 1, 1, :) = data_mat(size(data_mat, 1), :);
    
    % import to letswave
    [dataset(a).header, dataset(a).data, ~] = RLW_import_MAT_variable(data_mat_lw);

    % adjust the header
    dataset(a).header.name = data_name;
    dataset(a).header.xstep = 1/output_info(subject_idx).recording.SR;
    for c = 1:length(dataset(a).header.chanlocs)
        dataset(a).header.chanlocs(c).labels = output_info(subject_idx).recording.channels{c};
    end
end

% backup original dataset
dataset_raw = dataset;

% provide update
fprintf('Done. %d datasets imported.\n', a)
fprintf('\n')

% save info structure and move on
save(output_file, 'output_info', '-append');
cd(folder.data)
clear a b c data_name file2import data_mat data_mat_lw

%% pre-processing 1
% ----- section input -----
param.suffix = {'dc'};
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% loop through datasets
for d = 1:length(dataset)
    % provide update
    fprintf('processing dataset ''%s'':\n', dataset(d).header.name)

    % select the data
    lwdata.header = dataset(d).header;
    lwdata.data = dataset(d).data;

    % assign electrode coordinates
    fprintf('assigning electrode coordinates...')
    option = struct('filepath', sprintf('%s\\letswave 7\\res\\electrodes\\spherical_locations\\Standard-10-20-Cap81.locs', folder.toolbox), ...
        'suffix', '', 'is_save', 0);
    lwdata = FLW_electrode_location_assign.get_lwdata(lwdata, option);
    if d == 1
        output_info(subject_idx).preprocessing(1).process = sprintf('electrode coordinates assigned (standard 10-20-cap81)');
        output_info(subject_idx).preprocessing(1).date = sprintf('%s', date);
    end

    % remove DC + linear detrend
    fprintf('removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', param.suffix{1}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('DC correction + linear detrend');
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
    end
    
    % update the dataset
    dataset(d).header = lwdata.header;
    dataset(d).data = lwdata.data;
end

% provide update
fprintf('Done.\n')
fprintf('\n')

% save info structure and move on
save(output_file, 'output_info', '-append');
clear param d lwdata option 

%% ERPs: identify triggers
% ----- section input -----
param.datasets = [3, 4];
param.trigchan = {'force' 'force'};
param.trigname = {'pinprick' 'pinprick'};
param.threshold = [50 50];
param.isi = [10 10];
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% indentify and save triggers for all ERP datasets
for a = 1:length(param.datasets) 
    % provide update
    fprintf('processing dataset ''%s'':\n', dataset(param.datasets(a)).header.name)

    % identify trigger channel
    for b = 1:length(dataset(param.datasets(a)).header.chanlocs)
        if strcmp(param.trigchan{a}, dataset(param.datasets(a)).header.chanlocs(b).labels)
            trigchan = b;
        end
    end
    if ~exist('trigchan')
        fprintf('ERROR: no matching channel was found for chosen trigger source!')
    end

    % extract data, compute amplitude change and rectify the signal
    trigdata = squeeze(dataset(param.datasets(a)).data(1, trigchan, 1, 1, 1, :))';
    for c = 2:length(trigdata)
        trigdata_change(c-1) = trigdata(c) - trigdata(c-1);
    end
    trigdata_change = sqrt(trigdata_change.^2);

    % identify triggers
    trigpos = [];
    for d = 1:length(trigdata_change)
        if trigdata_change(d) > param.threshold(a) 
            if length(trigpos) == 0 
                trigpos(end+1) = d;
            elseif abs(trigpos(end)-d)*dataset(param.datasets(a)).header.xstep > param.isi(a)
                trigpos(end+1) = d;
            end
        end
    end

    % create the event file for letswave
    for e = 1:length(trigpos)
        dataset(param.datasets(a)).header.events(e).code = param.trigname{a};
        dataset(param.datasets(a)).header.events(e).latency = trigpos(e)*dataset(param.datasets(a)).header.xstep;
        dataset(param.datasets(a)).header.events(e).epoch = e;
    end

    % launch the figure
    x = [1:length(trigdata)]*dataset(param.datasets(a)).header.xstep;
    fig = figure(figure_counter);
    set(fig, 'units','normalized','outerposition',[0 0 1 1])
    xlim([x(1), x(end)]); ylim([-250, 750])
    xlabel('time (s)'); ylabel('amplitude (\muV)')
    title(sprintf('%s', dataset(param.datasets(a)).header.name))
    set(gca, 'fontsize', 14)
    hold on

    % plot the trigger data
    P(1) = plot(x, trigdata);
    P(2) = plot(x(1:end-1), trigdata_change);
    P(3) = plot(x, repelem(param.threshold(a), 1, length(x)));
    set(P(1), 'color', [0.0745    0.6235    1.0000], 'linewidth', 1.3)
    set(P(2), 'color', [0.0039    0.3216    0.5294], 'linewidth', 1.3)
    set(P(3), 'linewidth', 1.8)

    % plot the triggers
    for e = 1:length(trigpos)
        xline(trigpos(e)*dataset(param.datasets(a)).header.xstep, ...
            ':', 'color', [0.8000    0.0157    0.0157], 'linewidth', 2);
        hold on
    end
    
    % save figure
    save(fig, sprintf('%s\\figures\\triggers_%s_%d.png', folder.output, output_info(subject_idx).ID, a))
    
    % update figure counter
    figure_counter = figure_counter + 1;
    
    % update info structure
    if a == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: triggers identified based on channel data');
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
    end
end

% provide update
fprintf('Done.\n')
fprintf('\n')

% save info structure and move on
save(output_file, 'output_info', '-append');
clear param a b c d e trigchan trigdata trigdata_change trigpos x fig P

%% ERP: pre-processing 2
% ----- section input -----
param.datasets = [3, 4];
param.suffix = {'reref' 'bandpass' 'epoch', 'dc'};
param.ref = {'Cz','C3','C4','T7','T8','Pz'};
param.bandpass = [0.1 30];
param.epoch = [-0.5, 1];
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% loop through datasets
for d = 1:length(param.datasets)
    % select the data
    lwdata.header = dataset(param.datasets(d)).header;
    lwdata.data = dataset(param.datasets(d)).data;
    
    % provide update
    fprintf('processing dataset ''%s'':\n', lwdata.header.name)

    % re-reference to common average
    fprintf('re-referencing to common average...')
    option = struct('reference_list', {param.ref}, 'apply_list', {{lwdata.header.chanlocs(1:length(lwdata.header.chanlocs)-1).labels}}, ... 
        'suffix', param.suffix{1}, 'is_save', 0);
    lwdata = FLW_rereference.get_lwdata(lwdata, option);
    if d == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: re-referenced to common average: %s', [param.ref{:}]);
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
    end

%     % bandpass
%     fprintf('applying Butterworth bandpass filter...')
%     option = struct('filter_type', 'bandpass', 'high_cutoff', param.bandpass(2),'low_cutoff', param.bandpass(1),...
%         'filter_order', 4, 'suffix', param.suffix{2}, 'is_save', 0);
%     lwdata = FLW_butterworth_filter.get_lwdata(lwdata, option);
%     if d == 1
%         output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: bandpass filtered [%.1f %.1f]Hz - Butterworth, 4th order', param.bandpass(1), param.bandpass(2));
%         output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
%     end
    
    % segment to epochs
    fprintf('epoching from %d to %d ms relative to stimulus...\n', param.epoch(1)*1000, param.epoch(2)*1000)
    option = struct('event_labels', lwdata.header.events(1).code, 'x_start', param.epoch(1), 'x_end', param.epoch(2), ...
        'x_duration', param.epoch(2)-param.epoch(1), 'suffix', param.suffix{3}, 'is_save', 0);
    lwdata = FLW_segmentation.get_lwdata(lwdata, option);
    if d == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: segmented [%d %d]ms relative to stimulus', param.epoch(1)*1000, param.epoch(2)*1000);
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
    end
    
    % remove DC + linear detrend
    fprintf('removing DC and applying linear detrend...')
    option = struct('linear_detrend', 1, 'suffix', param.suffix{4}, 'is_save', 1);
    lwdata = FLW_dc_removal.get_lwdata(lwdata, option);
    if d == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: DC correction + linear detrend on single epochs');
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
    end

%     % notch + save
%     fprintf('applying FFT notch filter...')
%     option = struct('filter_type', 'notch', 'notch_fre', param.notch, 'notch_width', 2, 'slope_width', 2,...
%         'harmonic_num', 2,'suffix', param.suffix{3},'is_save', 1);
%     lwdata = FLW_FFT_filter.get_lwdata(lwdata, option);
%     if d == 1
%         output_info(subject_idx).preprocessing(end+1).process = sprintf('FFT notch filtered at 50 Hz');
%         output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
%     end

    % update the dataset
    dataset(param.datasets(d)).header = lwdata.header;
    dataset(param.datasets(d)).data = lwdata.data;
end

% provide update
fprintf('Done.\n')
fprintf('\n')

% save info structure and move on
save(output_file, 'output_info', '-append');
clear param d lwdata option 

%% visual check - eyeblink removal
% open letswave 6 for visual check
addpath(genpath([folder.toolbox '\letswave 6']));
% letswave

%% encode manual steps and finish processing
% ----- section input -----
param.datasets = [3, 4];
param.suffix = {'blinks' 'avg'};
% ------------------------- 
% add letswave 7 to the top of search path
addpath(genpath([folder.toolbox '\letswave 7']));

% loop through datasets
for d = 1:length(param.datasets)
    % select the data
    lwdata.header = dataset(param.datasets(d)).header;
    lwdata.data = dataset(param.datasets(d)).data;
    
    % provide update
    fprintf('processing dataset ''%s'':\n', lwdata.header.name)
    
    % get the info about removed epochs
    prompt = {'removed epochs:'};
    dlgtitle = sprintf('%s', dataset(param.datasets(d)).header.name);
    dims = [1 35];
    definput = {''};
    input = inputdlg(prompt,dlgtitle,dims,definput);
    clear prompt dlgtitle dims definput
    
    % fill in the structure
    if d == 1
        output_info(subject_idx).preprocessing(end+1).process = sprintf('ERP: epochs with blink artifact removed (visual inspection)');
        output_info(subject_idx).preprocessing(end).date = sprintf('%s', date);
        output_info(subject_idx).preprocessing(end).epochs{d, :} = num2str(input{1});
    else
        output_info(subject_idx).preprocessing(end).epochs{d, :} = num2str(input{1});
    end

    % average across trials
    fprintf('averaging across trials...')
    option = struct('operation', 'average', 'suffix', [param.suffix{2} ' ' param.suffix{1}], 'is_save', 1);
    lwdata = FLW_average_epochs.get_lwdata(lwdata, option);
end

% save info structure and move on to next subject
save(output_file, 'output_info', '-append');
folder.input = uigetdir(pwd, 'Coose the new input folder');             
clear d dataset side subject_idx input param lwdata option

%% helper lines
% load again the 'dataset' variable
file2import = dir(sprintf('*ep*S005*.mat'));
seq = [3, 4, 1, 2, 7, 8, 5, 6];
for d = 1:8
    load(file2import(seq(d)).name)
    dataset(d).data = data;
end
clear file2import seq d header data

% save again some data from 'dataset' variable
for d=3:length(dataset)
    LW_save(dataset(d).header.name, '', dataset(d).header, dataset(d).data)
end

% get rid of events with no-zero latency
for d = 1:length(dataset)
    % identify non-zero epochs
    for e = 1:length(dataset(d).header.events)
        if dataset(d).header.events(e).latency == 0
            idx(e) = false;
        else
            idx(e) = true;
        end
    end
    % remove these events 
    dataset(d).header.events(idx) = [];
    clear idx
    % save header to letswave
    header = dataset(d).header;
    save(sprintf('%s.lw6', dataset(d).header.name), 'header')
end
clear d e header

