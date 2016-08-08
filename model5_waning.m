%Model 5 - Waning model. Fixed: mu, gamma. Note compartment order is different for this model than the others!

% NIMBioS Cholera Project
% Created by: Brad Ochocki
% Modified by Marisa Eisenberg 1/14/2016

function [dxdt] = model5_waning(t, x, p)

% Number of S compartments
bins  = size(x,1) - 2;

% Define parameters
mu    = 1/(55*365);
gamma = 0.25;
alpha = p(4)*(bins-1);
ksi   = p(3);
B_I   = p(1);
B_W   = p(2);
u1    = p(5); % parameter for f(i)
u2    = p(6); % parameter for g(i)

% Equations for transitions due to infections (f(i)) and water (g(i))
f = exp(-linspace(0,10,bins)*u1);  % transitions to infections due to I
g = exp(-linspace(0,10,bins)*u2);  % transitions to infections due to W

% Set up x vector
W = x(1);
I = x(2);
S = x(3:end);

% Initiate dxdt
dxdt = nan(size(x,1),1);

% Define equations for I and W
dxdt(2) = sum(B_I*g*I.*S') + sum(B_W*f*W.*S') - gamma*I - mu*I;
dxdt(1) = ksi*(I-W);

% Define equation for S1
dxdt(3) = mu + alpha*S(2) - B_I*g(1)*S(1)*I - B_W*f(1)*S(1)*W - mu*S(1);

% Define equations S2:S(N-1)
for j = 2:(bins-1)
    dxdt(2+j) = alpha*S(j+1) - alpha*S(j) - B_I*g(j)*S(j)*I - B_W*f(j)*S(j)*W - mu*S(j);
end

% Define equation for SN
dxdt(end) = gamma*I - alpha*S(end) - B_I*g(end)*S(end)*I - B_W*f(end)*S(end)*W - mu*S(end);