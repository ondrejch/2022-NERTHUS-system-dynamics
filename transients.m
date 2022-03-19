%%% NERTHUS Scalable MSR System
%%% Authors - Nicholas Dunkle & Visura Pathirana 
%%% Project advisor - Dr. Ondrej Chvala

%% File Explanation
% This file collects transient information and all parameters for sim

% Transient types in this file:
% Reactivity insertions (benchmarks included)
% Load Following
% Pump trips & decay heat removal functionality

% Important variables to check before running simulation:
% fuel_type
% react_type
% load_follow_type

%% Basic Simulation Parameters
simtime = 10000;                                                           %Simulation time [s]
ts_max = 1e-1;                                                             %Maximum timestep [s] 
P=557;                                                                     %Operational thermal power [MW]
    % Primarily designed for 557MW. Works best scaled from ~100MW to ~1000MW

%% Fuel Types & Depletion
% fuel_type = 235; for FLibe with U235
% fuel_type = 233; for FLiBe with U233
% fuel_type = 123; for fuel with depletion accounting (requires depletion results file)
fuel_type = 123;                                                           
depletion_time = 0; %[days] integers only. Used for fuel type == 123

%% Reactivity Insertions

react_type = 0; % Set as 0 to use manual react data below. Set as 1-3 for benchmarks (next section)

%%% Source Step Reactivity Insertions & Sinusoidal Reactvity Insertions 
sourcedata = [0 0 0];                                                      %Neutron source insertions [abs]
sourcetime = [0 1000 2500];                                                %Neutron source insertion time [s]
source = timeseries(sourcedata,sourcetime);                                %Defining source timeseries  
% Reactivity step insertions
if react_type == 0
    reactdata = [0 0 0];                                                   %Reactivity insertions [abs]
    reacttime = [0 500 1000];                                              %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
% Reactivity sinusoidal insertions
omega          = 10.00000;                                                 %Frequncy of the sine wave [rad]
sin_mag        = 0;                                                        %Amplitude of the sine wave [abs]
dx             = round((2*pi/omega)/25, 2, 'significant');                 %Size of the time step [s]
    % To turn off sinusoidal insertion, set sin_mag = 0

%% Benchmarks
% Used if react_type == 1, 2, or 3

pcm2beta = 10^(-5);       % Converts pcm to absolute reactivity [unitless]
dollar_reactivity = 492;  % (pcm) Reactivity necessary for prompt critical insertion
bench1 = dollar_reactivity * 0.1 * pcm2beta;    % First benchmark insertion.    Prompt subcritical (1/10th PC)
bench2 = dollar_reactivity * pcm2beta;          % Second benchmark insertion.   Prompt critical (PC)
bench3 = dollar_reactivity * 2 * pcm2beta;      % Third benchmark insertion.    Prompt supercritical (2x PC)

if react_type == 1
    reactdata = [0 bench1];                                                %Reactivity insertions [abs]
    reacttime = [0 2000];                                                  %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
if react_type == 2
    reactdata = [0 bench2];                                                %Reactivity insertions [abs]
    reacttime = [0 2000];                                                  %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
if react_type == 3
    reactdata = [0 bench3];                                                %Reactivity insertions [abs]
    reacttime = [0 2000];                                                  %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end

%% Load Following

load_follow_type = 0; % Set as 0 to turn off. Set as 1 for load following

LF_rate = 10;                                      % Rate of power demand change [% nominal per minute]
LF_min = 5;                                        % Minimum power demand [% nominal power]
LF_min2 = LF_min/100;                              % Minimum power demand [power / nominal power]
LF_time = ((100 - LF_min)*60)/(LF_rate);           % Time of demand change [sec]
rise = linspace(LF_min2, 1, 200);                  % Rise magnitude
fall = linspace(1, LF_min2, 200);                  % Fall magnitude
LF_t1 = 2000 + LF_time;                            % End time of power demand fall
LF_t2 = LF_t1 + (2*LF_time);                       % Start time of power demand rise
LF_t3 = LF_t2 + LF_time;                           % End time of power demand rise
falltime1 = linspace(2000, LF_t1, 200);            % Fall duration
risetime1 = linspace(LF_t2, LF_t3, 200);           % Rise duration

% Load follow scenario feedwater flow rate demand
fdw_nom = P/2.8935; % The nominal full power feedwater flow rate
if load_follow_type == 0
    fdwdata = [fdw_nom fdw_nom];    
    fdwtime = [0 simtime];      
    fdw_demand = timeseries(fdwdata,fdwtime);
end
if load_follow_type == 1
    fdw_min = (557*0.05)/2.8935;                        % Minimum flow rate
    fdw_rise = linspace(fdw_min, fdw_nom, 200);         % Flow rate increase
    fdw_fall = linspace(fdw_nom, fdw_min, 200);         % Flow rate decrease
    fdwdata = [fdw_nom fdw_fall fdw_rise fdw_nom];      % Flow rates
    fdwtime = [0 falltime1 risetime1 simtime];          % Change times
    fdw_demand = timeseries(fdwdata,fdwtime);           % Feedwater demand
end

%% Pump Trips
    % To turn off pump trip, set trip time >>> simtime
Trip_P_pump = 20000000;                                                    %Time at which primary pump is tripped [s]
Trip_S_pump = 20000000;                                                    %Time at which secondary pump is tripped [s]
Trip_T_pump = 20000000;                                                    %Time at which secondary pump is tripped [s]
Trip_UHX    = 20000000;                                                    %Time at which ultimate heat exchanger will be cut off [s]
    % To turn off UHX trip, set trip time >>> simtime 
    % Reminder: UHX must be attached for this function to work
    % Default is OTSG, not UHX

%% Decay Heat Removal System (DHRS)

%%% DHRS_MODE = 1; a sigmoid based DHRS (Normal DHRS)
%%% DHRS_MODE = 2; a square pulse based DHRS (Broken DHRS)
%%% DHRS_MODE = 1; allows modifications to sigmoid behavior using parameters with Normal DHRS in parameter file
%%% DHRS_MODE = 2; allows cold slug insertions
DHRS_MODE = 1; 
DHRS_time=20000000;                                                        %Time at which DRACS will be activated [s]
    % To turn off DRACS, set DHRS_time >>> simtime

%%% Only for DHRS_MODE = 1
DHRS_Power= P*(0.10);                                                      %Maximum power that can be removed by DHRS
Power_Bleed= P*(0.00);                                                     %Some power can be removed from DRACS even when its not activated 
    % Make sure Power_Bleed = 0 if you don't want it

%%% Only for DHRS_MODE = 2
SlugDeltaTemp = 30;                                                        %Temperature drop by broken DHRS [deg. C]
Slug_duration = 10;                                                        %Duration of slug [s]

%% Parameters 

run('parameters.m');    % Runs the parameters file to get all the relavent simulation variables into the workspace
