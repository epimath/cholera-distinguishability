%Model 2 - Gamma model. Fixed: mu, gamma.

% NIMBioS Cholera Project
% Created by: Mike Kelly
% Modified by Marisa Eisenberg 1/14/2016

function dxdt = model2_gamma(t, x, params)

number = length(x(:,1))-4;

mho = 1/(55*365);
bI = params(1);
bW = params(2);
gamm = 0.25;
xi = params(3);
alph = params(4)*(number+1);
k = params(5);

S = x(1);
I = x(2);
W = x(3);
R = x(4);
R_I = zeros(number,1);

for j=1:number
    R_I(j,1) = x(j+4);
end

dxdt = zeros(4+number,1);
dxdt(1)=mho-bI*S*I-bW*S*W-mho*S+alph*R_I(number);
dxdt(2)=bI*I*S+bW*S*W-gamm*I-mho*I;
dxdt(3)=xi*(I-W);
dxdt(4)=gamm*I-mho*R-alph*R;

for j=1:number
    if j==1
    dxdt(4+j)=alph*R-mho*R_I(j)-alph*R_I(j);
    else
    dxdt(4+j)=alph*R_I(j-1)-mho*R_I(j)-alph*R_I(j);
    end
end
