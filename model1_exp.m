%Model 1 - Exponential model. Fixed: mu, gamma

% NIMBioS Cholera Project
% Simulation Model 1: SIWR
% Created by: Karen Hamre
% Last udpated: August 19, 2013
% Modified by Marisa Eisenberg 1/14/2016

function [dxdt]=model1_exp(t,x,params)
mu=1/(55*365);              %params(1), can set to 0 when assume constant birth and death rate in population; 
betaI=params(1);            %params(2), direct transmission;
betaW=params(2);            %params(3), water-borne transmission;
gamma=0.25;                 %params(4), recovery rate;
xi=params(3);               %params(5), rate of decay;
alpha=params(4);            %params(6), loss of immunity;
k = params(5);              %params(7), measurement error;

%Use when simulating data for endemic settings with seasonal transmission, comment out for epidemics
%betaWseas = betaW/2*(1+ sin(2*t*pi/365));  


S=x(1);
I=x(2);
W=x(3);
R=x(4);


dxdt=zeros(4,1);

% Comment out one below depending on transmission status
dxdt(1)=mu-betaI*S*I-betaW*S*W-mu*S+alpha*R;  %For epidemics
%dxdt(1)=mu-betaI*S*I-betaWseas*S*W-mu*S+alpha*R;  %For endemic seasonal transmission

% Comment out one below depending on transmission status
dxdt(2)=betaI*S*I+betaW*S*W-gamma*I-mu*I;     %For epidemics
%dxdt(2)=betaI*S*I+betaWseas*S*W-gamma*I-mu*I;     %For endemic seasonal transmission

dxdt(3)=xi*(I-W);
dxdt(4)=gamma*I-mu*R-alpha*R;


