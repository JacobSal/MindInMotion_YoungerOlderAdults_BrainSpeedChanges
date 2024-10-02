function X = cpu_par_func(data, numWorkers)

    Nsig = size(data,1);
    Nsamp = size(data,2);

    fb = cwtfilterbank('SignalLength',Nsamp,'Voices',10, 'Wavelet', 'amor',...
        'SamplingPeriod', hours(1/60), 'PeriodLimits', [hours(1) hours(32)]);
    Ns = length(fb.Scales);
    X = zeros([Ns,Nsamp,Nsig],'single');
    data = single(data);
    parfor (ii = 1:Nsig, numWorkers)
        X(:,:,ii) = abs(fb.wt(data(ii,:)));
    end
    
end

