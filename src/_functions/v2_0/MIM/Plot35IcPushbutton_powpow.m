function Plot35IcPushbutton_powpow(EEG,numIC,IC_IDs)
% hObject    handle to Plot35IcPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% EEG = evalin('base', 'EEG');
freqs = EEG.etc.PowPowCAT.freqs;

figure
set(gcf, 'color', [0.66 0.76 1])
numPanels = min([size(EEG.etc.PowPowCAT.covMatrix,3) numIC]);
for panelIdx = 1:numPanels
    
    % Plot correlation matrix
    if     numPanels <= 6
        subplot(2,3,panelIdx)
    elseif numPanels <= 12
        subplot(3,4,panelIdx)
    elseif numPanels <= 24
        subplot(4,6,panelIdx)
    else
        subplot(5,7,panelIdx)
    end

    imagesc(EEG.etc.PowPowCAT.covMatrix(:,:,panelIdx), [-0.8 0.8]);
    colormap(jet);
    axis xy
    axis square
    
    % Zoom in a little bit
    gcaPosition = get(gca, 'position');
    gcaPosition([3 4]) = gcaPosition([3 4])*1.1;
    set(gca, 'position', gcaPosition);
    tickLabels = round(freqs(10:10:length(freqs))*10)/10;
    tickLabels(tickLabels>10) = round(tickLabels(tickLabels>10));
    set(gca, 'XTick', [], 'XTickLabel', [],...
        'YTick', 10:10:length(freqs), 'YTickLabel', tickLabels, 'fontsize', 8)
    title(['IC ' num2str(IC_IDs(panelIdx))], 'fontsize', 12)
    if panelIdx == size(EEG.icaweights,1)
        return
    end
end
