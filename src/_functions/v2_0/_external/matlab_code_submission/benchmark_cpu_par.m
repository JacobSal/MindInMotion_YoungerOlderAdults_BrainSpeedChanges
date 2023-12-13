%% Setup
myCluster = parcluster('local');
finalNumWorkers = 64;
parpool(finalNumWorkers);

numDays = 10*30;
minsInDay = 60*24;

t = 1:1:minsInDay*numDays;
full_data = zeros(1000,length(t));

% Want 1 cycle every 1440 mins/24 hours
y = sin((2.0*pi*t)/minsInDay);
noiseAmplitude = .1;
for i= 1:size(full_data,1)
    full_data(i,:) = y + noiseAmplitude * randn(1, length(y));
end

wanted_num_samples = [[1:99] [100:25:1000]];

all_data = {};
for i = 1:length(wanted_num_samples)
    all_data{i} = full_data(1:wanted_num_samples(i),:);
end

timed_arr = zeros(length(wanted_num_samples),10);

for i = 1:length(wanted_num_samples)
    disp(i);
    for j = 1:5
        tic;
        X = cpu_par_func(all_data{i},finalNumWorkers);
        timed_arr(i,j)= toc;
        X = [];
        writematrix(timed_arr,'cpu_par_benchmarks.csv');
    end
end
delete(gcp('nocreate'));