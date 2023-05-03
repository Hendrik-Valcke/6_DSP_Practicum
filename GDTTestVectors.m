time = 0:10:(15*60);

intervals = length(time);



start_load = 0;
max_load = 13;
load = [];
for i = start_load:max_load % Generate data like 0 0 .. 0 10 10 .. 10 20 20 .. 20
    xi = repelem(i*10,6); 
    load = cat(2, load, xi);   
end

% Match time vector
if length(load) < length(time)
    load(end+1:length(time)) = 130;
end

% define parameters
n_values = intervals; % number of pulse rate values to generate
min_value = 60; % minimum pulse rate value
max_value = 150; % maximum pulse rate value
max_step_size = 5; % maximum step size between adjacent values

% generate random step sizes
step_sizes = randi([-max_step_size, max_step_size], 1, n_values-1);

% generate pulse rate values with increasing or decreasing step sizes
pulse_rate = [min_value cumsum(step_sizes) + min_value];

% plot the generated pulse rate values
plot(time,pulse_rate);
hold on
plot(time,load)


