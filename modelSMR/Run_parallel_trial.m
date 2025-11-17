clear; clc; close all;

%% === Number of Monte Carlo simulations ===
N = 14;   % you can change this to 10, 100, etc.

% Random seed (optional, for reproducibility)
rng(42);

%% === Define random sampling ranges for each input ===
ranges = struct( ...
    'LOCA_time',   [1 1200], ...
    'LOOP_time',   [1 1200], ...
    'LHS_time',    [1 1200], ...
    'MSLB_time',   [1 1200], ...
    'MSLBH2_time', [1 1200], ...
    'LOOPH2_time', [1 1200], ...
    'ACS_time',    [1 1200], ...
    'ACS_rtime',   [30   90], ...
    'EDG_time',    [1 1200], ...
    'EDG_rtime',   [10   180], ...
    'pdp_time',    [1  1200], ...
    'pdp_rtime',   [10   180] ...
);

%% === Randomly sample input values within defined ranges ===
LOCA_time   = randRange(ranges.LOCA_time,   N);
LOOP_time   = randRange(ranges.LOOP_time,   N);
LHS_time    = randRange(ranges.LHS_time,    N);
MSLB_time   = randRange(ranges.MSLB_time,   N);
MSLBH2_time = randRange(ranges.MSLBH2_time, N);
LOOPH2_time = randRange(ranges.LOOPH2_time, N);
ACS_time    = randRange(ranges.ACS_time,    N);
ACS_rtime   = randRange(ranges.ACS_rtime,   N);
EDG_time    = randRange(ranges.EDG_time,    N);
EDG_rtime   = randRange(ranges.EDG_rtime,   N);
pdp_time    = randRange(ranges.pdp_time,    N);
pdp_rtime   = randRange(ranges.pdp_rtime,   N);

%% === Run the parallel simulation ===
[T_W1, T_W2, T_W3, T_W4] = Simulation_model_hydrogen_parallel( ...
    LOCA_time, LOOP_time, LHS_time, ...
    MSLB_time, MSLBH2_time, LOOPH2_time, ...
    ACS_time, ACS_rtime, ...
    EDG_time, EDG_rtime, ...
    pdp_time, pdp_rtime);

%% === Plot results ===
%t = 1:size(T_W1,2); % time vector (assuming tsim steps)

%figure('Name','Hydrogen Parallel Simulation Results','Color','w');
%subplot(2,2,1); plot(t, T_W1'); title('T\_W1'); xlabel('Time'); ylabel('Temperature');
%subplot(2,2,2); plot(t, T_W2'); title('T\_W2'); xlabel('Time'); ylabel('Temperature');
%subplot(2,2,3); plot(t, T_W3'); title('T\_W3'); xlabel('Time'); ylabel('Temperature');
%subplot(2,2,4); plot(t, T_W4'); title('T\_W4'); xlabel('Time'); ylabel('Temperature');

%sgtitle(sprintf('Hydrogen Parallel Simulation Results (%d runs)', N));

%% === Helper function ===
% function vals = randRange(range, N)
%     % Uniform random sampling in [min, max]
%     vals = range(1) + (range(2)-range(1)) * rand(1, N);
% end
