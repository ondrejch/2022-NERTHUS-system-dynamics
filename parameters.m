%%% NERTHUS Scalable MSR System
%%% Authors - Nicholas Dunkle & Visura Pathirana 
%%% Project advisor - Dr. Ondrej Chvala

%% Simulation Parameter File Explanation
% This script calculates and initializes simulation parameters
% User inputs for transient designs are in transients.m

% WARNING: Do not mess with this file unless you know what you are doing.
% The file designed to be modified for most use cases is transients.m

% For more information on the model, seek the following journal articles:
%%% Modeling Methodology: "Molten Salt Reactor Power Plant System Dynamic Modeling using Nodal Methodology"
%%% NERTHUS Overview: "NERTHUS Thermal Molten Salt Reactor Dynamic Model in MATLAB-Simulink"
%%% Load Following with NERTHUS: "Effect of Xenon Removal Rate on Load Following in High Power Thermal Spectrum Molten-Salt Reactors (MSRs)"

%% Scaling Factors (SF)

%%% Core and Heat Exchanger Scaling Factor Calculations
Power_ratio_CHE = P/557; % The reactor power divided by what the power the components were initially designed for
SF_fdw_CORE = ((P-557)/557)*4209;   % Scaling factor for the fluid flow rate through the reactor core and primary loop
SF_fdw_PHX = ((P-557)/557)*1543;    % Scaling factor for the fluid flow rate in the secondary loop
SF_fdw_SHX = ((P-557)/557)*3379;    % Scaling factor for the fluid flow rate in the tertoary loop
SF_linear_CHE = Power_ratio_CHE;                % The scaling factor used for linear scaling (masses and volumes)
SF_length_CHE = (Power_ratio_CHE)^(1/3);        % The scaling factor used for lengths
SF_area_CORE  = (Power_ratio_CHE)^(2/3);        % The scaling factor used for areas
    % CHE is an acronym for "Core & Heat Exchangers"

%%% Steam Generator Scaling Factor Calculations
Power_ratio = P/750; % The reactor power divided by what the power the steam generator was initially designed for
SF_fdw = P/2.8935;   % Steam Generator feedwater flow rate
SF_linear = Power_ratio;                        % The scaling factor used for linear scaling (masses and volumes)
SF_length = (Power_ratio)^(1/3);                % The scaling factor used for lengths
SF_area = (Power_ratio)^(2/3);                  % The scaling factor used for areas

%% Circulation Parameters (Salt resident times) {tau_x*vdot_x=volume_x}
    % Tau is the resident time of a unit of fluid in that component
% Primary loop
tau_core = 15;              % [s]
tau_phx = 11.6;             % [s]
tau_DHRS = 1;               % [s]
tau_core_DHRS = 1;          % [s]
tau_DHRS_phx = 1;           % [s]
tau_phx_core = 1;           % [s]
tau_loop = tau_phx + tau_core_DHRS + tau_DHRS_phx + tau_phx_core;

% Secondary loop
tau_coolant_phx = 8;        % [s]
tau_coolant_shx = 14;       % [s]
tau_coolant_phx_shx = 1;    % [s]
tau_coolant_shx_phx = 1;    % [s]

% Tertiary loop
tau_HITEC_shx = 11;         % [s]
tau_HITEC_UHX = 10;         % [s]
tau_HITEC_shx_UHX = 1;      % [s]
tau_HITEC_UHX_shx = 1;      % [s]

vdot_f_fc_pc = 0.05; % Free convection flowrate percentage in the primary loop and secondary loop
K_pump = 10000000; % Pump coast down constant for both pumps (modeled as exp decay) 

% Fluid flow rates (kg/s)
mdot_fuel = 4209 + SF_fdw_CORE;   %Fuel salt mass flow rate (kg/s)    - Primary Loop
mdot_coolant = 1543 + SF_fdw_PHX; %Coolant salt mass flow rate (kg/s) - Secondary Loop
mdot_HITEC = 3379 + SF_fdw_SHX;   %Solar salt mass flow rate (kg/s)   - Tertiary Loop

%% Material Properties

%Material Densities (kg/m^3)
rho_fuel     = 2.14647E+03; % (partially enriched U-235)ORNL-TM-0728 p.8 2.243E+03; % (Th-U) density of fuel salt (kg/m^3) ORNL-TM-0728 p.8
rho_graphite = 1.860E3;     % graphite density (kg/m^3) ORNL-3812 p.77, ORNL-TM-0728 p.87
rho_coolant  = 1.922e3;     % coolant salt density (kg/m^3) ORNL-TM-0728 p.8
rho_tube     = 8.7745E+03;  % (kg/m^3) density of INOR-8 ORNL-TM-0728 p.20
rho_HITEC    = 1680;        % (kg/m^3) density of Hitec at 550C

%Material Specific Heat Capacities (MJ/kg-C)
scp_fuel     = 1.9665E-3;   % specific heat capacity of fuel salt (MJ/kg-C) ORNL-TM-0728 p.8
scp_graphite = 1.773E-3;    % cp_g/m_g; % graphite specific heat capacity (MW-s/kg-C) ORNL-TM-1647 p.3
scp_coolant  = 8.83E-3;     % specific heat capacity of coolant (MJ/(kg-C) ORNL-TM-0728 p.8
scp_tube     = 5.778E-04;   % specific heat capacity of tubes (MJ/(kg-C)) ORNL-TM-0728 p.20  
scp_HITEC    = 0.00154912;  % specific heat capacity of HITEC salt (MJ/(kg-C))

%% Modified Point Kinetics Module (mPKE)

n_frac0 = 1;             % initial fractional nuetron density n/n0 [n/cm^3/s]

%%%U-235 data
if fuel_type == 235
    Lam  = 2.400E-04;    % Mean generation time {ORNL-TM-1070 p.15 for U235}
    lam  = [1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00];
    beta = [0.000223,0.001457,0.001307,0.002628,0.000766,0.00023]; % U235
    a_f  = -8.71E-05;    % U235 (drho/°C) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
    a_g  = -6.66E-05;    % U235 (drho/°C) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101
end

%%%U-233 data
if fuel_type == 233
    Lam    = 4.0E-04;    % mean generation time ORNL-TM-1070 p.15 U233
    lam    = [1.260E-02, 3.370E-02, 1.390E-01, 3.250E-01, 1.130E+00, 2.500E+00]; %U233
    beta   = [0.00023,0.00079,0.00067,0.00073,0.00013,0.00009]; % U233
    a_f    = -11.034E-5; % U233 (drho/°C) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
    a_g    = -05.814E-5; % U233  (drho/°C) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101
end

%%%U-235 data with depletion
if fuel_type == 123
    run('read_depletion.m');
    depl_index = find(range_depl  == depletion_time);
    Lam    = depl_Lam(depl_index);
    lam    = [1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00];
    beta   = depl_beta(depl_index,:);
    a_f    = 1E-5*depl_fuel_temp_coef(depl_index);
    a_g    = 1E-5*depl_grap_temp_coef(depl_index);
end    

beta_t = sum(beta); % Total delayed neutron fraction

% reactivity change in going from stationary to circulating fuel
C0(1)  = ((beta(1))/Lam)*(1.0/(lam(1) - (exp(-lam(1)*tau_loop) - 1.0)/tau_core));
C0(2)  = ((beta(2))/Lam)*(1.0/(lam(2) - (exp(-lam(2)*tau_loop) - 1.0)/tau_core));
C0(3)  = ((beta(3))/Lam)*(1.0/(lam(3) - (exp(-lam(3)*tau_loop) - 1.0)/tau_core));
C0(4)  = ((beta(4))/Lam)*(1.0/(lam(4) - (exp(-lam(4)*tau_loop) - 1.0)/tau_core));
C0(5)  = ((beta(5))/Lam)*(1.0/(lam(5) - (exp(-lam(5)*tau_loop) - 1.0)/tau_core));
C0(6)  = ((beta(6))/Lam)*(1.0/(lam(6) - (exp(-lam(6)*tau_loop) - 1.0)/tau_core));

% Reactivity loss from circulating fuel 
bterm    = 0;
for i = 1:6
    bterm = bterm + beta(i)/(1.0 + ((1.0-exp(-lam(i)*tau_loop))/(lam(i)*tau_core)));
end
rho_0    = beta_t - bterm;      % reactivity loss due to delayed neutrons born outside of the core

%% Benchmarks
% Used if react_type == 1, 2, or 3 (transients.m file)

dollar = (beta_t - rho_0) * 1e5;     % 1$ of Reactivity expressed in pcm. Modified by circulatory loss of DNPs
pcm2beta = 10^(-5); % Converts pcm to absolute reactivity [unitless]
bench1 = dollar * 0.1 * pcm2beta;    % First benchmark insertion.    Prompt subcritical (1/10th PC)
bench2 = dollar * pcm2beta;          % Second benchmark insertion.   Prompt critical (PC)
bench3 = dollar * 2 * pcm2beta;      % Third benchmark insertion.    Prompt supercritical (2x PC)
bench4 = 100 * pcm2beta;             % Fourth benchmark 

if react_type == 1
    reactdata = [0 bench1];                                                %Reactivity insertions [abs]
    reacttime = [0 benchmark_start];                                       %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
if react_type == 2
    reactdata = [0 bench2];                                                %Reactivity insertions [abs]
    reacttime = [0 benchmark_start];                                       %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
if react_type == 3
    reactdata = [0 bench3];                                                %Reactivity insertions [abs]
    reacttime = [0 benchmark_start];                                       %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end
if react_type == 4
    reactdata = [0 bench4];                                                %Reactivity insertions [abs]
    reacttime = [0 benchmark_start];                                       %Reactivity insertion time [s]
    react = timeseries(reactdata,reacttime);                               %Defining source timeseries
end

%% Core Design Parameters
W_f      = mdot_fuel;           % Fluid flow rate (kg/s)
m_f      = mdot_fuel*tau_core;  % fuel mass in core (kg)
nn_f     = 2;                   % number of fuel nodes in core model
mn_f     = (m_f/nn_f);          % fuel mass per node (kg)
f_g_Ratio = 0.03;               % Fuel to graphite mass ratio

% Core Upflow
m_g      = m_f/f_g_Ratio;       % graphite mass (kg)

% mcp is material mass * heat capacity of the material
mcp_g1   = m_g*scp_graphite;    % mcp of graphite per lump (MW-s/C)
mcp_f1   = mn_f*scp_fuel;       % mcp of fuel salt per lump (MW-s/C)
mcp_f2   = mn_f*scp_fuel;       % mcp of fuel salt per lump (MW-s/C)

% hA is the heat transfer coefficient * heat transfer area
hA_fg    = 12.5325*SF_area_CORE; % HTC x heat transfer area (MW/C) calculated from MSRE scaling 

k_g      = 0.07;     % fraction of total power generated in the graphite  ORNL-TM-0728 p.9
k_1      = 0.5;      % fraction of heat transferred from graphite which goes to first fuel lump
k_2      = 0.5;      % fraction of heat transferred from graphite which goes to second fuel lump
k_f      = 0.93;     % fraction of heat generated in fuel - that generated in external loop ORNL-TM-0728 p.9
k_f1     = k_f/nn_f; % fraction of total power generated in lump f1
k_f2     = k_f/nn_f; % fraction of total power generated in lump f2

Tf_in  = 616.6;      % Core inlet  fuel temp (C)
T0_f2  = 681.8;      % Core outlet fuel temp (C)

T0_f1  = Tf_in + (T0_f2-Tf_in)/2;
T0_g1  = T0_f1 + (k_g*P/hA_fg); 

%% Primary Heat Exchanger (PHX) parameters

Tc_in  = 586; % Core Inlet  Temp (C)
Tc_out = 627; % Core Outlet Temp (C)

% UA is calculated from log mean temperature drop as follows,
delta_T1 = T0_f2 - Tc_out;
delta_T2 = Tf_in - Tc_in;
delta_T_lm = (delta_T1 - delta_T2)/log(delta_T1/delta_T2);

UA_phx = 5*P/delta_T_lm;

Wp_phx1   =  mdot_fuel;                 % primary fluid mass flow rate in PHX, 1/6 of primary flow rate, table 2 ORNL-TM-3832 (kg/s)
Ws_phx1   =  mdot_coolant;              % secondary fluid mass flow rate in PHX, table 4 ORNL-TM-3832 (kg/s)

mp_phx = mdot_fuel*tau_phx;             % Total mass of fuel in primary side of PHX
ms_phx = mdot_coolant*tau_coolant_phx;  % Total mass of coolant in secondary side of PHX

mp_phx_nm = mp_phx/4;                   % Individual primary node mass of fuel in PHX
ms_phx_nm = ms_phx/4;                   % Individual secondary node mass of coolant in PHX

mt_phx1_node = 1323.6;%*SF_linear_PHX;
mp_phx1   = [mp_phx_nm, mp_phx_nm, mp_phx_nm, mp_phx_nm];   % PHX1 primary node masses (kg)
mt_phx1   = [mt_phx1_node, mt_phx1_node];                   % PHX1 tube node masses (kg)
ms_phx1   = [ms_phx_nm, ms_phx_nm, ms_phx_nm, ms_phx_nm];   % PHX1 secondary node masses (kg)

% mcp = (mass)*(specific heat capacity)
mcp_p_phx1 = scp_fuel*mp_phx1;          % PHX1 primary nodes
mcp_t_phx1 = scp_tube*mt_phx1;          % PHX1 tube nodes
mcp_s_phx1 = scp_coolant*ms_phx1;       % PHX1 secondary nodes

UA_phx_nm = UA_phx/4; %Individual node hA from log mean temperature calculation 

% hA = (heat transfer coefficient)*(area)
hA_p_phx1 = [UA_phx_nm,UA_phx_nm,UA_phx_nm,UA_phx_nm];      % for primary/tube interface
hA_s_phx1 = [UA_phx_nm,UA_phx_nm,UA_phx_nm,UA_phx_nm];      % for tube/secondary interface

%Initial temperature
T0_p_phx1 = [663.1, 644.3, 630.3, 616.6];   % PHX1 primary nodes
T0_t_phx1 = [639.3, 612.7];                 % PHX1 tube nodes
T0_s_phx1 = [595.1, 603.8, 615.6, 627.3];   % PHX1 secondary nodes

%% Secondary Heat Exchanger (SHX) parameters

Tss_in = 507.3;  % PHX Secondary Side Inlet  Temp (C)
Tss_out = 613.7; % PHX Secondary Side Outlet Temp (C)

% UA is calculated from log mean temperature drop as follows,
delta_T1_shx = Tc_out - Tss_in;
delta_T2_shx = Tss_out - Tc_in;
delta_T_lm = (delta_T1_shx - delta_T2_shx)/log(delta_T1_shx/delta_T2_shx);

UA_shx = 40*P/delta_T_lm;

Wp_shx1 = mdot_coolant; % kg/s
Ws_shx1 = mdot_HITEC;   % kg/s

mp_shx = mdot_coolant*tau_coolant_shx;  % Total mass of coolanr in primary side of SHX
ms_shx = mdot_HITEC*tau_HITEC_shx;      % Total mass of HITEC in secondary side of SHX

mp_shx_nm = mp_shx/4;   % per node
ms_shx_nm = ms_shx/4;   % per node

mt_shx1_tube = 1940;% * SF_linear_SHX;
mp_shx1   = [mp_shx_nm, mp_shx_nm, mp_shx_nm, mp_shx_nm];   % SHX1 primary nodes
mt_shx1   = [mt_shx1_tube, mt_shx1_tube];                   % SHX1 tube nodes
ms_shx1   = [ms_shx_nm, ms_shx_nm, ms_shx_nm, ms_shx_nm];   % SHX1 secondary nodes

% mcp = (mass)*(specific heat capacity)
mcp_p_shx1 = scp_coolant*mp_shx1;   % SHX1 primary nodes
mcp_t_shx1 = scp_tube*mt_shx1;      % SHX1 tube nodes
mcp_s_shx1 = scp_HITEC*ms_shx1;     % SHX1 secondary nodes

UA_shx_nm = UA_shx/4; %Individual node hA from log mean temperature calculation 

hA_p_shx1 = [UA_shx_nm,UA_shx_nm,UA_shx_nm,UA_shx_nm];
hA_s_shx1 = [UA_shx_nm,UA_shx_nm,UA_shx_nm,UA_shx_nm];

% Initial temperatures (C)
T0_p_shx1 = [641.8, 656.3, 621.4, 586.4];   % SHX1 primary nodes
T0_t_shx1 = [646.6, 609.8];                 % SHX1 tube nodes
T0_s_shx1 = [598.2, 689.1, 651.4, 613.7];   % SHX1 secondary nodes

%% Pipe nodes for decay heat
mn_piDHRS = mdot_fuel*tau_DHRS;
m_f_t = mdot_fuel*(tau_core + tau_phx + tau_DHRS); % Total mass of fuelsalt in primary loop

%% Decay heat data
% this splits the fission products into 
Fiss_factor = 2.464783802008740e-03; % (rel to power) heat per fission relative to power rating [calculated]

% fission yields for each of the three groups lumped with heating factor
gamma0  = 9.635981959409105e-01;
gamma1  = 3.560674858154914e-02;
gamma2  = 7.950554775404400e-04;

% decay constants for each of the gorups
lambda0 = 0.0945298;
lambda1 = 0.00441957;
lambda2 = 8.60979e-05;

%% DRACS Parameters

% Normal DHRS
% Power_Bleed= P*(0.01); %Some power will removed from DRACS even when its not used 
Epsilon=1E-3;
DHRS_TIME_K=10;

% Broken DHRS
rm_power = W_f*scp_fuel*SlugDeltaTemp;
slug_end = DHRS_time+Slug_duration;
slug_heat_rm = [0 rm_power 0 0 ];
slugtime = [0 DHRS_time slug_end simtime];
slugger = timeseries(slug_heat_rm,slugtime);

%% Steam Generator (SG)

%%% Heat capacities and densities
Cp_p        = 0.00154912;  % Hitec specific heat [[MJ/(kg/C)]
Cp_w        = 456.056e-6;  % [[MJ/(kg/C)]  % inconel                        
rho_p       = rho_HITEC;   % [kg/m^3]      % Hitec salt          
rho_w       = 8425.712;    % [kg/m^3]      % inconel                        

%%% Steam Generator Tube Side

N           = 6546*SF_linear;      % 6546 Number of Tubes     
L_ft        = 28*SF_length;        % [ft] need to use in a few places
L           = 8.5344*SF_length;    % [m]             % Active Tube Length
L_b         = 2.3645*SF_length;    % [m]             % Boiling Length    
L_s         = 4.723*SF_length;     % [m]             % Steam Length    
L_sc        = 1.4469*SF_length;    % [m]             % Subcooled Length  
D_ot        = 0.015875;            % [m]             % outer tube diameter            [ft]
T_th        = 0.0008636;           % [m]             % tube thickness                 [ft]
D_it        = 0.0141478;           % [m]             % internal tube diameter         [ft]
R_it        = 0.0070739;           % [m]             % internal tube radius           [ft] 
R_ot        = 0.0079375;           % [m]             % outer tube radius              [ft]
A_sit       = 2483.0632*SF_area;   % [m^2]           % total surface area inner tube  [ft^2]
A_sot       = 2786.2020*SF_area;   % [m^2]           % total surface area outer tube  [ft^2]

%%% Initial State Calculations
M_stm   = 0.018;       % [kg/mol]        18.0000       ; % Molar weight of steam          [lbm/lb-mol]
P_table = 11.0:0.1:12.5; %Pressure;
Ts_avg = (459 + 411)/2 + 273;
T_table = [];
Hfg_table = [];
hs_table = [];
for PPP = 11.0:0.1:12.5
[T_sat,hf,hg,kf,kg] = hsat(PPP);
[~, hss, ~, ~] = hsh(PPP, Ts_avg);
T_table = [T_table, T_sat];
hs_table = [hs_table, hss];
Hfg_table = [Hfg_table, hg-hf];
end
T_table = (T_table - 273); % Saturated temperature;
Hfg_table = Hfg_table*1e-6;
hs_table = hs_table*1e-6;
a = polyfit(P_table, T_table, 1);
X_5 = a(2); K_5 = a(1);
b = polyfit(P_table, Hfg_table, 1);
X_4 = b(2); K_4 = b(1);
c = polyfit(P_table, hs_table, 1);
dHs_dPs = c(1);

%%%%%%%%%% Temperature %%%%%%%%%%
T_p1        = 605.5;     % [C]           % Primary Coolant Temperature node
T_p2        = 579.6;     % [C]           % Primary Coolant Temperature node
T_p3        = 557.5;     % [C]           % Primary Coolant Temperature node
T_p4        = 537.3;     % [C]           % Primary Coolant Temperature node
T_p5        = 523.9;     % [C]           % Primary Coolant Temperature node
T_p6        = 507.3;     % [C]           % Primary Coolant Temperature node
T_w1        = 598.4;     % [C]           % Temperature for wall node 1     
T_w2        = 557.3;     % [C]           % Temperature for wall node 2     
T_w3        = 414.7;     % [C]           % Temperature for wall node 3     
T_w4        = 407.0;     % [C]           % Temperature for wall node 4     
T_w5        = 424.3;     % [C]           % Temperature for wall node 5     
T_w6        = 382.7;     % [C]           % Temperature for wall node 6     
T_s1        = 594.7;     % [C]           % Temperature for superheated node
T_s2        = 537.6;     % [C]           % Temperature for superheated node
T_sc2       = 250.1;     % [C]           % Temperature for subcooled node2 
T_fw        = 180.0;     % [C]           % Feedwater temperature           
T_sat       = 327.8;     % [C]           % Saturation temperature of the s[F]
T_pin       = 613.7;     % [C]           % primary inlet temperature      [F]

%%%%%%%%%% Pressure %%%%%%%%%%
P_p1        = 0;
P_p2        = 0;
P_p3        = 0;
P_p4        = 0;
P_p5        = 0;
P_p6        = 0;
P_sat       = 12.5; % [MPa]     % Saturation Pressure       
P_set       = 12.5; % [MPa]     % Steam Pressure Setpoint     
P_ss0       = 12.5; % [MPa]     % Pressure superheated steam lump
P_s         = 12.5; % [MPa]     % Pressure superheated steam lump
P_sc        = 12.5; % [MPa]     % Pressure subcooled initial     

%%%%%%%%%% Mass Flow Rate %%%%%%%%%%
W_3 = mdot_HITEC; % HITEC Salt flow rate in the tertiary loop
W_p0        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate inlet
W_p1        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate node 1
W_p2        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate node 2
W_p3        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate node 3
W_p4        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate node 4
W_p5        = W_3; % 3779.937;    % [kg/s]        % Mass Flow Rate outlet
W_fw        = SF_fdw;             % [kg/s]        % Mass flow rate of feedwater

%%%%%%%%%% Area %%%%%%%%%%
A_pw1       = 186.2298*SF_area;    % [m^2]         % heat transfer area of Node 1
A_pw2       = 186.2298*SF_area;    % [m^2]         % heat transfer area of Node 2
A_pw3       = 869.0721*SF_area;    % [m^2]         % heat transfer area of Node 3
A_pw4       = 869.0721*SF_area;    % [m^2]         % heat transfer area of Node 4
A_pw5       = 186.2298*SF_area;    % [m^2]         % heat transfer area of Node 5
A_pw6       = 186.2298*SF_area;    % [m^2]         % heat transfer area of Node 6
A_ws1       = 208.9652*SF_area;    % [m^2]         % Area of heat transfer node 1
A_ws2       = 208.9652*SF_area;    % [m^2]         % Area of heat transfer node 2
A_ws3       = 975.1707*SF_area;    % [m^2]         % Area of heat transfer node 3
A_ws4       = 975.1707*SF_area;    % [m^2]         % Area of heat transfer node 4
A_ws5       = 208.9652*SF_area;    % [m^2]         % Area of heat transfer node 5
A_ws6       = 208.9652*SF_area;    % [m^2]         % Area of heat transfer node 6
    % These last three are not scaled because the tube diameters are not scaled
A_s         = 1.2245;      % [m^2]          % Cross sectional secondary flow
A_p         = 1.0034;      % [m^2]          % Primary flow area  
A_w         = 0.2666;      % [m^2]          % Cross section for the tube

K_b         = 12.834854076292217;  % 
K_1         = 11.632248097704855;  %
K_sc        = -17.615718797133468; % 
dHs_dPs     = -0.042433622114357;  %

Z_ss        = 0.76634;      % Compressibility factor at 570 K, 60 atm
R           = 8.314462E-6;  % [MJ/mol-�C] % Universal gas constant 

%%%%%%%%%% HEAT TRANSFER COEFFICIENTS %%%%%%%%%%
    % HTCs are not scaled because the proportionality between the heat flux and thermo dynamic driving force is constant
h_s_shx     = 9101.343578935E-6;
h_pw        = h_s_shx; %9101.343578935E-6; % [MW/m^2-C]      1.8070       ; % Effective primary to wall heat [BTU/s-ft^2-F]
h_ws        = 5732.378387768E-6;    % [MW/m^2-C]      0.6672       ; % htc wall to steam node         [BTU/s-ft^2-F]
h_wb        = 13349.334646671E-6;   % [MW/m^2-C]      2.16647      ; % htc wall to boil node          [BTU/s-ft^2-F]
h_wsc       = 8385.005556375E-6;    % [MW/m^2-C]      1.18         ; %0.147 htc wall to subcooled node     [BTU/s-ft^2-F]

%%%%%%%%%% Steam Enthalpy Calculation%%%%%%%%%%
temp_table = linspace(180, 650, 500);   % Temp range in C
pres_inp = linspace(120, 160, 500);     % Pressure range in Bar
pres_table = pres_inp;
m = length(temp_table);
n = length(pres_inp);
H_s_table = zeros(m, n); % matrix of zeros to store result
for i=1:m
    for j=1:n
        H_s_table(i, j) = XSteam('h_pT', pres_inp(j), temp_table(i))/1e3; % in MJ/kg
    end
end
H_s_table2 = zeros(1,m);
for k=1:m
    H_s_table2(k) = XSteam('h_pT', 125, temp_table(k))/1e3;
end
%%%%%%%%%% Specific Heat Calculation%%%%%%%%%%
temp_table2 = linspace(325, 700, 500); % Temp range in C
pres_table2 = linspace(120, 160, 500); % Pressure range in Bar
m2 = length(temp_table2);
n2 = length(pres_table2);
Cp_table = zeros(m2, n2); % matrix of zeros to store result
for i2=1:m2
    for j2=1:n2
        Cp_table(i2, j2) = XSteam('Cp_pT', pres_table2(j2), temp_table2(i2))/1e3; % in MJ/kg
    end
end
%    Feedwater    %%%%%%%%%%%%%%%%%%%%
Cal_temp = linspace(180,630,450);
Cal_Cp1 = zeros(1,400);
Cal_temp1 = linspace(180, 200,400);
region_fdw = linspace(T_fw, 200, 400);
for z = 1:400
    Cal_Cp1(z) = XSteam('Cp_pT', 125, Cal_temp1(z));
end
Cp_1 = mean(interp1(Cal_temp1, Cal_Cp1, region_fdw))/1e3;
%    Subcooled    %%%%%%%%%%%%%%%%%%%%
Cal_Cp2 = zeros(1,400);
region_sc = linspace(200, 327.8, 400);
Cal_temp2 = linspace(200, 327.8, 400);
for z = 1:400
    Cal_Cp2(z) = XSteam('Cp_pT', 125, Cal_temp2(z));
end
Cp_2 = mean(interp1(Cal_temp2, Cal_Cp2, region_sc))/1e3;
%      Steam      %%%%%%%%%%%%%%%%%%%%
Cal_Cp3 = zeros(1,400);
region_s2 = linspace(327.82, 535, 400);
Cal_temp3 = linspace(327.82, 535, 400);
for z = 1:400
    Cal_Cp3(z) = XSteam('Cp_pT', 125, Cal_temp3(z));
end
Cp_3 = mean(interp1(Cal_temp3, Cal_Cp3, region_s2))/1e3;
%    Region 4
Cal_Cp4 = zeros(1,400);
region_s1 = linspace(535, 605, 400);
Cal_temp4 = linspace(535, 605, 400);
for z = 1:400
    Cal_Cp4(z) = XSteam('Cp_pT', 125, Cal_temp4(z));
end
Cp_4 = mean(interp1(Cal_temp4, Cal_Cp4, region_s1))/1e3;
% Cp_4 = Cp_3;

%%%%%%%%%% Density Calculation %%%%%%%%%%
%       Feedwater    %%%%%%%%%%%%%%%%%%%%
Cal_rho1 = zeros(1,400);
for z = 1:400
    Cal_rho1(z) = XSteam('rho_pT', 125, Cal_temp1(z));
end
rho_fw = mean(interp1(Cal_temp1, Cal_rho1, region_fdw));
%rho_fw = XSteam('rho_pT', 125, 180);
%       Subcool    %%%%%%%%%%%%%%%%%%%%
Cal_rho2 = zeros(1,400);
for z = 1:400
    Cal_rho2(z) = XSteam('rho_pT', 125, Cal_temp2(z));
end
rho_sc = mean(interp1(Cal_temp2, Cal_rho2, region_sc));
rho_b = (XSteam('rho_pT', 125, 327.81) + XSteam('rho_pT', 125, 327.82))/2;
%       Steam    %%%%%%%%%%%%%%%%%%%%
Cal_rho3 = zeros(1,400);
Cal_temp_steam = linspace(328, 595, 400);
region_steam = linspace(328, 595, 400);
for z = 1:400
    Cal_rho3(z) = XSteam('rho_pT', 125, Cal_temp_steam(z));
end
rho_s = mean(interp1(Cal_temp_steam, Cal_rho3, region_steam));

M_initial   = L_b*rho_b*A_s + L_s*rho_s*A_s + L_sc*A_s*rho_sc;

%% Xenon reactivity effects

% fission yield
gamma_I      = 0.06390;     % fission yield for I-135 and short lived parents (Tc & Sb) for U-235 [www-nds.iaea.org]
gamma_Xe     = 0.0025761;   % fission yield for Xe-135 for U-235 [www-nds.iaea.org]
gamma_Pm     = 0.010816;    % fission yield for Pm-149 and short lived parents (Nb & Pr) for U-235 [www-nds.iaea.org]
gamma_Sm     = 0.0;         % fission yield for Sm-149 for U-235 [www-nds.iaea.org]

% decay constants
I_135_lam    = log(2)/(6.58*60*60);  % (s^-1) radioactive decay constant for I-135 [www-nds.iaea.org]
Xe_135_lam   = log(2)/(9.14*60*60);  % (s^-1) radioactive decay constant for Xe-135 [www-nds.iaea.org]
Pm_149_lam   = log(2)/(53.08*60*60); % (s^-1) radioactive decay constant for Pm-149 [www-nds.iaea.org]

% absorption cross sections
sig_Xe       = 2.66449e-18;  % (cm^2) mircoscopic cross-section for Xe (n,gamma) reaction [www-nds.iaea.org]
sig_Sm       = 4.032e-20;    % (cm^2) mircoscopic cross-section for Sm (n,x) reaction [www-nds.iaea.org]

% xenon removal rate 
lam_bubble = 0.00500; % [s^-1] 0.005 is moderate; 0.02039 (MSRE's rate) is high
%%% For low removal rates, saturation of the xenon concentration is
%%% necessary to have the model successfully initialized. Controls for the
%%% saturation are lam_b_upp, lam_b_down, and saturation_end.

% Fuel parameters
mol_wt_fuel  = 64.64; % (g) weight of a mol of fuel salt
mol_den_fuel = 0.0519; % (mol/cm3) molar denisty of fuel salt
U_den        = 1.5632e20; % (#U/cm^3)
U_sig        = 5.851e-22; % (cm^2) Microscopic fission cross section of the core

% Neutronics parameters
%%% Determinging average Sig_f and Sig_a  by weighting flux 
Sig_f = 1.37520e-4 * (9.48937e+21 / (9.48937e+21 + 1.05375e+22)) + 1.50804e-3 * (1.05375e+22 / (9.48937e+21 + 1.05375e+22));
Sig_a = 1.48040e-3 * (9.48937e+21 / (9.48937e+21 + 1.05375e+22)) + 2.59018e-3 * (1.05375e+22 / (9.48937e+21 + 1.05375e+22));
phi_0 = 3.1849E+14; % From Serpent: total integrated flux divided by volume of lattice region

% initial conditions
I_0          = gamma_I*Sig_f*phi_0/I_135_lam; % (#/cm^3) steady state iodine concetration
Xe_0         = ((gamma_Xe+gamma_I)*Sig_f*phi_0)/(Xe_135_lam + sig_Xe*phi_0 + lam_bubble); % (#/cm^3) steady state xenon concetration
Pm_0         = gamma_Pm*Sig_f*phi_0/Pm_149_lam; % (#/cm^3) steady state Promethium concetration
Sm_0         = ((gamma_Pm+gamma_Sm)*Sig_f)/(sig_Sm);

% Xenon saturation block controls 
lam_b_upp = Xe_0*(1.01);    % Upper Xe concentration limit
lam_b_down = Xe_0*(0.99);   % Lower Xe concentration limit
saturation_end = 11000;     % [sec] time where saturation releases

% reactivity from Fission product poisons
rho_Xe       = (Xe_0*sig_Xe)/Sig_a;
rho_Sm       = (Sm_0*sig_Sm)/Sig_a;


%% Alternative Ultimate Heat Sink
% UHX have three nodes that runs HITEC and remove demand power

n_UHX = 3;

mn_UHX = (tau_HITEC_UHX*mdot_HITEC)/n_UHX;

tau_ostg_shx = tau_HITEC_UHX_shx;

