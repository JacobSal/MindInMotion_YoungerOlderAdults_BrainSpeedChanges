% Fooof test
pyversion
eeglab

settings = struct();  % Use defaults
f_range = [1, 40];

%% Load data
load('M:\liu.chang1\STUDY-preprocess-HY_20220831\BATCH-3-Epoch\ERSP_Plots\Cluster_7\readSPEC_1_Terrain.mat');

% Input should be in linear spacing
specdata_nolog = 10.^(specdata{1}/10);
% Run FOOOF
return_model = true;
fooof_results = fooof(specfreqs, specdata_nolog(:,1), f_range, settings, return_model);

% Print out the FOOOF Results
fooof_results

% Notes:
% # Initialize a FOOOF model object with defined settings
% fm = FOOOF(peak_width_limits=[1.0, 8.0], max_n_peaks=6, min_peak_height=0.1,
%            peak_threshold=2.0, aperiodic_mode='fixed'
% Use settings:
%   settings        = struct, can optionally include:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
fm = py.fooof.FOOOF(settings.peak_width_limits, ...
                        settings.max_n_peaks, ...
                        settings.min_peak_height, ...
                        settings.peak_threshold, ...
                        settings.aperiodic_mode, ...
                        settings.verbose);
fm.fit(py.numpy.array(specfreqs), py.numpy.array(specdata_nolog),f_range);

%% plot the result
log_freq = 0;
single_fig = 1;

fooof_plot(fooof_results,log_freq,single_fig);
set(gcf,'color','white');
xlabel('Frequency(Hz)');
ylabel('Power');

% Second plot, remove the initial
init_flat_spec =  fooof_results.power_spectrum - fooof_results.ap_fit;
figure();
plot(fooof_results.freqs, init_flat_spec,'color','black');
% Find peak 
set(gcf,'color','white');
xlabel('Frequency(Hz)');
ylabel('Power');

%% Test fooof group analysis
settings.peak_width_limits = [1.0, 8.0];
settings.max_n_peaks = 8;
settings = fooof_check_settings(settings);

% fg group doesn't work in matlab. Has to loop through 
for i = 1:size(specdata_nolog,2)
    fooof_group_results{i} = fooof(specfreqs, specdata_nolog(:,i), f_range, settings, return_model);
end

% Output from the fooof results
% Aperiodic_params =  [Offset, (Knee), Exponent]
% Peak_params = [CF, PW, BW]; Center freq, Power, Bandwidth
% Gaussian fits : each row is a gaussian as mean, height, standard
% deviation

% Define frequency bands of interest
%% Plot distribution of aperiodic params (exp), central frequency, and goodness of fit
for i = 1:size(fooof_group_results,2)
    fooof_group_results_org(i).aperiodic_exp = fooof_group_results{i}.aperiodic_params(2);
    fooof_group_results_org(i).central_freq  = fooof_group_results{i}.peak_params(:,1);
    fooof_group_results_org(i).r_squared     = fooof_group_results{i}.r_squared;
end
figure();set(gcf,'color','white');
subplot(2,2,1)
plot(ones(size(fooof_group_results,2),1),[fooof_group_results_org(:).aperiodic_exp],'.');
ylabel('Aperodic exponent');

subplot(2,2,2)
histogram(vertcat(fooof_group_results_org(:).central_freq));
xlabel('Central frequency');ylabel('# Peaks')

subplot(2,2,3)
plot(ones(size(fooof_group_results,2),1),[fooof_group_results_org(:).r_squared],'.');
ylabel('Central frequency');

%{
fg = py.fooof.FOOOFGroup(settings.peak_width_limits,settings.max_n_peaks)

% Fit FOOOF model across the matrix of power spectra
f_range = py.list(f_range);
fg_results = fg.fit(py.numpy.array(specfreqs), py.numpy.array(specdata_nolog'),f_range)

% Create and save out a report summarizing the results across the group of power spectra
fg.save_report()

% Save out FOOOF results for further analysis later
fg.save('fooof_group_results', True)
%}