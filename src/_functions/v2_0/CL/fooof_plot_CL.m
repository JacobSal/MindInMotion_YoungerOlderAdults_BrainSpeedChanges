% fooof_plot() - Plot a FOOOF model.
%
% Usage:
%   >> fooof_plot(fooof_results)
%
% Inputs:
%   fooof_results   = struct of fooof results
%                       Note: must contain FOOOF model, not just results
%   log_freqs       = boolean, whether to log frequency axis
%                       Note: this argument is optional, defaults to false
%

function fooof_plot_CL(fooof_results, log_freqs, single_fig, legend_on, mode)
    % if single_fig, creates a new figure
    %   assign false if using subplots
    if ~exist('single_fig', 'var')
        single_fig = true;
    end
    %% Data Checking

    if ~isfield(fooof_results, 'freqs')
       error('FOOOF results struct does not contain model output.')
    end

    %% Set Up

    if exist('log_freqs', 'var') && log_freqs
        plt_freqs = log10(fooof_results.freqs);
    else
        plt_freqs = fooof_results.freqs;
    end

    % Plot settings
    lw = 2.5;
    
    theta_band = [4, 8];
    alpha_band = [8 12];
    beta_band  = [12 30];

    %% Create the plots
    if single_fig
%         figure()
    end
    hold on

    if strcmp(mode, 'flatten_curve')
        data = plot(plt_freqs, fooof_results.power_spectrum - fooof_results.ap_fit , 'black');hold on;
        xline(4,'color','b');xline(8,'color','b');xline(12,'color','b');xline(30,'color','b');
    else
        % Plot the original data
        data = plot(plt_freqs, fooof_results.power_spectrum, 'black');

        % Plot the full model fit
        model = plot(plt_freqs, fooof_results.fooofed_spectrum, 'red');

        % Plot the aperiodic fit
        ap_fit = plot(plt_freqs, fooof_results.ap_fit, 'b--');

    end

    %% Added by Chang
    % Add the central frequency and power cursor
    switch mode
        case 'line'
            for i = 1:size(fooof_results.peak_params,1)
                xline(fooof_results.peak_params(i,1),'g','linewidth',2);
            end
        case 'outline'
            for i = 1:size(fooof_results.peak_params,1)
                peak = fooof_results.gaussian_params;
                peak_range = [peak(i,1) - peak(i,3)*3, peak(i,1) + peak(i,3)*3];
                peak_line = fooof_results.ap_fit + peak(i,2)*exp(-(fooof_results.freqs - peak(i,1)).^2/(2*peak(i,3).^2));
                peak_freqs = fooof_results.freqs(fooof_results.freqs > peak_range(1) & fooof_results.freqs < peak_range(2));
                peak_line = peak_line(fooof_results.freqs > peak_range(1) & fooof_results.freqs < peak_range(2));
                
                plot(peak_freqs,peak_line,'color','gree','linewidth',2);
            end
        case 'dot'
            for i = 1:size(fooof_results.peak_params,1)
                ap_point = interp1(fooof_results.freqs,fooof_results.ap_fit,fooof_results.peak_params(i,1));
                plot([fooof_results.peak_params(i,1) fooof_results.peak_params(i,1)],[ap_point ap_point + fooof_results.peak_params(i,2)],...
                    'color','green','linewidth',2);
                plot([fooof_results.peak_params(i,1)],[ap_point + fooof_results.peak_params(i,2)],...
                    '.','color','green','linewidth',2);
            end
        case 'shade'
        case 'outline_dot'
            for i = 1:size(fooof_results.peak_params,1)
                peak = fooof_results.gaussian_params;
                peak_range = [peak(i,1) - peak(i,3)*3, peak(i,1) + peak(i,3)*3];
                peak_line = fooof_results.ap_fit + peak(i,2)*exp(-(fooof_results.freqs - peak(i,1)).^2/(2*peak(i,3).^2));
                peak_freqs = fooof_results.freqs(fooof_results.freqs > peak_range(1) & fooof_results.freqs < peak_range(2));
                peak_line = peak_line(fooof_results.freqs > peak_range(1) & fooof_results.freqs < peak_range(2));
                ap_point = interp1(fooof_results.freqs,fooof_results.ap_fit,fooof_results.peak_params(i,1));
                plot([fooof_results.peak_params(i,1) fooof_results.peak_params(i,1)],[ap_point ap_point + fooof_results.peak_params(i,2)],...
                    'color','green','linewidth',1);
                plot([fooof_results.peak_params(i,1)],[ap_point + fooof_results.peak_params(i,2)],...
                    '.','color','green','linewidth',2);

                plot(peak_freqs,peak_line,'color','gree','linewidth',2);
            end
    end
    %% Plot Settings

    % Apply general plot settings
    if ~strcmp(mode, 'flatten_curve')
        for plt = [data, model, ap_fit]
            set(plt, 'LineWidth', lw);

            % Set alpha value for model - in a wonky way, because Matlab
            %   Note: the '4' is magical and mysterious. No idea.

            model.Color(4) = 0.5;

            grid on
            if legend_on
                legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit')
            end
            if single_fig
                hold off
            end
        end
    else
        set(data, 'LineWidth', lw);
        grid on, hold on;
    end
end