%Model 3 - Asymptomatic model. Fixed: mu, gamma, q; beta_IA constrained.

% NIMBioS Cholera Project
% Created by: Elizabeth Lee
% Modified by Marisa Eisenberg 1/14/2016

function dxdt = model3_asymp_restrict (t, x, params)

%% Define starting parameters defined in main code
% parameter order for simdata_main, fit_main, and forecast_main .m files
mu = 1/(55*365); 
beta_IS = params(1)*(5/2);  %to give 5/8 baseline value
beta_IA = beta_IS(1)/4;
% note the above makes 
% 0.2*betaIS + 0.8*beta_IA
% = (1/5) betaIS + (4/5) betaIA
% = params(1)/2 + params(1)/2 = params(1), which is our identifiable combo
beta_W = params(2);
alpha_S = params(3);
alpha_A = params(4);
q = 0.2;
gamma = 0.25;
ksi = params(5);
k = params(6);

%% State variables
S = x(1); % susceptible
I_S = x(2); % symptomatically infected
I_A = x(3); % asymptomatically infected
R_S = x(4); % recovered from symptomatic infection
R_A = x(5); % recovered from asymptomatic infection
W = x(6); % concentration of pathogen in the water

%% Differential equations - epidemic
dxdt = zeros(6,1); % S, I_S, I_A, R_S, R_A, W
dxdt(1) = mu - beta_IS*S*I_S - beta_IA*S*I_A - beta_W*S*W - mu*S + alpha_S*R_S + alpha_A*R_A;
dxdt(2) = q*beta_IS*S*I_S + q*beta_IA*S*I_A + q*beta_W*S*W - gamma*I_S - mu*I_S;
dxdt(3) = (1-q)*beta_IS*S*I_S + (1-q)*beta_IA*S*I_A + (1-q)*beta_W*S*W - gamma*I_A - mu*I_A;
dxdt(4) = gamma*I_S - alpha_S*R_S - mu*R_S; 
dxdt(5) = gamma*I_A - alpha_A*R_A - mu*R_A;
dxdt(6) = ksi*(I_A + I_S - W);

%% Additional comments

% Parameter definitions:
% mu - rate of change for birth and death rates; set at zero for now
% beta_IS - transmission parameter for infection due to interaction with symptomatics
% beta_IA - transmission parameter for infection due to interaction with asymptomatics
% beta_W - transmission paratmeter for infection due to contact with water
% alpha_S - loss of immunity for symptomatic individuals
% alpha_A - loss of immunity for asymptomatic individuals
% q - fraction of new infections that are symptomatic (20% commonly used)
% gamma - rate of recovery among infected individuals
% ksi - combined- shedding of pathogen into water and rate of
% decay of pathogen already in water
% k - 

% King et al. 2008 states that estimates of the ratio of asymptomatic to symptomatic
% infections has ranged from 3 to 100, but their study suggests that this
% ratio is high.

% 1) mu = rate of change in population via birth and death rates
% 2-4) beta_IS, beta_IA, beta_W = transmission parameter for infections due
% to contact with individuals with symptomatic and
% asymptomatic infections and rate at which susceptibles become infected
% due to contact with pathogen in the water
% 5-6) alpha_S, alpha_A = loss of immunity for symptomatic and asymptomatic
% individuals
% 7) q = fraction of new infections that are symptomatic
% 8) gamma = rate of recovery of infecteds 
% 9) xi = combined parameter -- shedding of pathogen into water and rate of
% decay of pathogen already in water

% Note: beta_IS*S*IS is the change of susceptibles to infected due to
% interaction with symptomatic individuals


