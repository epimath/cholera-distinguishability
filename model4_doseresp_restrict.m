%Model 4 - Dose Response model. Fixed: mu, gamma, n, K.

% NIMBioS cholera project
% ODE file for model 4
% Created by: Segun Akinwumi
% Last modified: Sep 20, 2014
% Adapted Sep 2015 by MCE
% Modified by Marisa Eisenberg 1/14/2016


function [dy]=model4_doseresp_restrict(t,x,params)
mu=1/(55*365);
gamma=0.25;
n=1;
K=0.4;
betaI=params(1);
betaW=params(2);
xi=params(3);
alpha=params(4);
S=x(1);
I=x(2);
W=x(3);
R=x(4);
dy=zeros(4,1);
dy(1)=mu-betaI*S*I-betaW*S*W^n/(K+W^n)-mu*S+alpha*R;
dy(2)=betaI*S*I+betaW*S*W^n/(K+W^n)-gamma*I-mu*I;
dy(3)=xi*(I-W);
dy(4)=gamma*I-mu*R-alpha*R;

