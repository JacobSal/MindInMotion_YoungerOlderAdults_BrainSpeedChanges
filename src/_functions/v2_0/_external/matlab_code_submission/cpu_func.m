function X = cpu_func(data)

    Nsig = size(data,1);
    Nsamp = size(data,2);

    fb = cwtfilterbank('SignalLength', Nsamp, 'Voices', 10, 'Wavelet', 'amor',...
        'SamplingPeriod', hours(1/60), 'PeriodLimits', [hours(1) hours(32)]);
    X = zeros([Ns,size(data,2),Nsig],'single');
    data = single(data);

    for ii = 1:Nsig
        cfs = abs(fb.wt(data(ii,:)));
        X(:,:,ii) = abs(cfs);
    end
end

