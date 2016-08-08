% Main Code - simulate data & fit with all models
% 1/14/2016 Marisa Eisenberg (marisae@umich.edu)

clear

%% Notes

% Files this needs: all the ode functions listed in models, as well as:
% cost_ML.m, fisher.m (if you want FIMs), AngolaData.m (if using data from
% Angola epidemic).

% How to use this code: 
%    - All the main settings you need are in the Model and Data Setup
%      section below. For most things, you can edit that section and then
%      just run this code, assuming you have the files listed above.
%    - However, for plotting the forecasting and Angola data results, there
%      are extra code files (labeled "mini") which do the plots in that
%      case. (This is silly and should really just be integrated in here or
%      turned into it's own plot&save script, but I'm lazy so here we are.)
%    - Be sure to actually check all the settings in the Model and Data
%      Setup section! This section is super long, but we did a lot of 
%      different variations (e.g. combinations of noise/cost function/starting
%      parameters/etc.), so it's not necessarily the case that changing one
%      thing to say 'poisson noise' means the corresponding cost function
%      or other settings also change.
%    - Currently, the Model and Data Setup section is set to simulate data
%      with normal (Gaussian) noise and estimate all the models from that 
%      using 'naive' starting parameters. 



%% Model and Data Setup

%%%% Misc Naming Stuff for Saving Files

% Name for saving files
filename = 'SimData_normnoise_naive_080716';

% Param plots & latextable
%   This is just labels for the sections at the end that keep track of the results for plotting and table-making.
paramnames = {'beta_i';'beta_w';'xi';'alpha';'k'};
latexparamnames = {'$\beta_i$';'$\beta_w$';'$\xi$';'$\alpha$';'$k$'};
extracols = [0 1 1 0]; %pois norm epidemic informed
noisestart = {'normal','naive'}; %'none','poisson','normal'; 'informed,'naive'



%%%% Data Setup

% Angola data
%   If using Angola data, comment/skip the data simulation section--and for 
%   plotting and saving, use mini_AngolaPlotSave.m
% AngolaData; 
% datasets = {data};

% Simulated data fitting
tspan = 0:1:100; %100 is duration of simulated data we've been using
% tspan = 0:1:161; %161 is duration of Angola data
% tspan = 0:1:3*365; %3 year - multi-season data to see multiple waves

% Forecasting
% tspan = 0:1:10;
% tspan = 0:1:30;
% tspan = 0:1:50;


% Data Noise Type
% noise = @(y) y;                 % no noise
% noise = @(y) poissrnd(y);       % poisson noise
noise = @(y) normrnd(y,0.1*y);  % normal noise with sigma = 0.1*y



%%%% Fitting Setup
% Measurement equation
yfcn = @(p,x) x(:,2)/p; 

% Cost function for fitting
% costfcn = @(data,y) (y)'*ones(length(y),1) - data'*log(y); % Poisson ML
% costfcn = @(data,y) (y - data)'*(y - data); %OLS 
costfcn = @(data,y) ((y - data)./sqrt(data))'*((y - data)./sqrt(data)); %WLS weight sigma^2 = data 
% costfcn = @(data,y) ((y - data)./(0.1*data))'*((y - data)./(0.1*data)); %WLS weight sigma^2 = (0.1*data)^2

% AIC
% aicfcn = @(data,p,cost) 2*length(p) + 2*cost + 2*sum(log(factorial(round(data)))); %Poisson AIC
% aicfcn = @(data,p,cost) 2*length(p) + 2*cost; %Truncated Poisson AIC - dropping the giant factorial term, because holy cats. We could pass the log into the factorial, but we really don't need this term anyhow, so ditching it.
aicfcn = @(data,p,cost) 2*length(p) + length(data)*log(2*pi()) + sum(log(data))+cost; %WLS weight sigma^2 = data
% aicfcn = @(data,p,cost) 2*length(p) + length(data)*log(2*pi()) + 2*sum(log(0.1*data))+cost; %WLS weight sigma = 0.1*data



%%%% Models Setup
% List of models
models = {@model1_exp,@model2_gamma,@model3_asymp_restrict,@model4_doseresp_restrict,@model5_waning};
modelnames = {'Exponential','Gamma','Asymptomatic','Dose Response','Waning'};



%%% Baseline parameters
% complete list of all parameters to be used
%           betaI,betaW, xi, alpha, k,      alpha_S,   alpha_A   u1,u2
allparams = [0.25 0.50 0.01 (1/365) 1/50000 1/(2*365) 1/(0.5*365) 1 1];
% use these to make the parameter sets for each model (same order as models)
baseparams = {allparams(1:5),allparams(1:5),...
    [allparams(1:2) allparams(6:7) allparams(3) 0.2*allparams(5)],... %note that order and k for asymp is different!
    allparams(1:5),[allparams(1:4) allparams(8:9) allparams(5)]};

dataparams = baseparams;

%%%% Naive starting values (within ±20% of true values):
baseparams = cellfun(@(p) rand(size(p)).*p*0.4 + 0.8*p,baseparams,'UniformOutput', false);



%%% Initial conditions
% Build functions that will generate ICs on the fly for each model 
%    num = total # of compartments for the model including SIWR 
genICs = @(data,k,num) [1-data(1)*k; data(1)*k; zeros(num-2,1)];
asympICs = @(data,k,num) [1-5*data(1)*k; data(1)*k; 4*data(1)*k; zeros(num-3,1)]; % S IS IA RS RA W
waningICs = @(data,k,num) [0; data(1)*k; 1-data(1)*k; zeros(num-3,1)];  %W I S0 S1..Sn
% Array them, same order as models. These are the same as the functions above, 
% but we're eliminating the num, so we don't have to enter it every time.
x0fcns = {@(data,k)genICs(data,k,4),@(data,k)genICs(data,k,12),...
            @(data,k)asympICs(data,k,6),@(data,k)genICs(data,k,4),@(data,k)waningICs(data,k,12)};

        
        
%%% R0 functions
basicR0 = @(params) (params(1) + params(2))/(0.25 + 1/(55*365));
doserespR0 = @(params) (params(1) + params(2)/0.4)/(0.25 + 1/(55*365));
r0fcns = {basicR0,basicR0,basicR0,doserespR0,basicR0};

% So, from here onwards, we'll use models, baseparams, x0fcns, costfcn, and yfcn for
% everything, and we can just loop through those as needed.

%% Simulate data for fitting

baseICdata = 0.01*50000;
datasets = {}; % yes, I should preallocate and it would be easy to. But meh. 
               % (Also get used to it, because it's going to happen a lot and a lot more egregiously if you keep reading...:P)

% figure(1); hold on

for i=1:length(models)
    [tsim,xsim] = ode45(models{i},tspan,x0fcns{i}(baseICdata,dataparams{i}(end)),[],dataparams{i});
    datasets = [datasets,noise(yfcn(dataparams{i}(end),xsim))];
%     plot(tsim,datasets{i},'o')
end

% datasets = datasets(5);

%% Fit to each data set with each model

paramests = {};
gofs = [];
flags = [];

for i=1:length(datasets)
    for j=1:length(models)
        [i j]
        [currparamests, currgof, success] = fminsearch(@(param) cost_ML(tspan, datasets{i}, models{j}, x0fcns{j}, param, yfcn, costfcn),baseparams{j},optimset('MaxFunEvals',8000,'MaxIter',8000));%,'TolX', 1.0000e-0,'TolFun',1.0000e-0));  % 'Display','iter',
        paramests{i,j} = abs(currparamests);
        gofs(i,j) = currgof;
        flags(i,j) = success;
    end
end




%%%%Comment from here downward if you are doing forecasting or Angola fits!
%%%% I should really either merge the forecasting/Angola stuff or make this
%%%% a separate function for everything. Thing to do for later!





%% Plot all fits
% note we can also integrate this with the fitting loops themselves, but
% doing it separately in case we want to plot from previous runs with new formatting, etc.

fits = {};

for i=1:length(datasets)
    figure(i)
    set(gca,'LineWidth',1,'FontSize',20,'FontName','Arial')
    hold on
    plot(tspan,datasets{i},'k.','LineWidth',2.5,'MarkerSize',12)
    for j=1:length(models)
        [tsim,xsim] = ode45(models{j},tspan,x0fcns{j}(datasets{i},paramests{i,j}(end)),[],paramests{i,j});
        fits{i,j} = yfcn(paramests{i,j}(end),xsim);
        plot(tsim,fits{i,j},'LineWidth',2.5)
    end
    title(strcat(modelnames{i},' Data'))
    ylabel('Cases');
    xlabel('Days');
    legend(['Data' modelnames]);
    
    %optional save figs
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    saveas(fig,strcat(filename,'_',modelnames{i},'Data'),'pdf')
    saveas(fig,strcat(filename,'_',modelnames{i},'Data'),'fig')
end

%% Save data


%%%% Matlab Data
% Workspace is pretty small, so just save everything in case we need it later.
save(strcat(filename,'.mat'))




%% Save Data for parameter plots & latex

%%%% Write data for parameter plots
paramplotdata = {};
fitplotdata = {};

psizer = ones(length(dataparams{1}),1);
dsizer = cell(size(tspan'));
msizer = cell(size(tspan'));
nssizer = {};
for i=1:length(tspan)
    nssizer = [nssizer; noisestart];
end

for i=1:length(datasets)
    for j = 1:length(models)
        if j == 3
            miniparams = [paramests{i,j}(1:2) paramests{i,j}(5) paramests{i,j}(3)*0.2+paramests{i,j}(4)*0.8 paramests{i,j}(end)/0.2];
        else
            miniparams = [paramests{i,j}(1:4) paramests{i,j}(end)];
        end  
        paramplotdata = [paramplotdata; num2cell([i*psizer extracols(1)*psizer extracols(2)*psizer j*psizer extracols(3)*psizer extracols(4)*psizer]) paramnames num2cell([dataparams{1}' miniparams'])];
        
        dsizer(:) = modelnames(i);
        msizer(:) = modelnames(j);
        fitplotdata = [fitplotdata; dsizer msizer nssizer num2cell(tspan') num2cell(fits{i,j})];
    end
    
    % add true data to fitplot data
    dsizer(:) = modelnames(i);
    msizer(:) = {'True Data'};
    fitplotdata = [fitplotdata; dsizer msizer nssizer num2cell(tspan') num2cell(datasets{i})];
end


%    write the little labels at the end - could merge this with above, but meh
for i=1:size(paramplotdata,1)
    paramplotdata{i,10} = strcat(num2str(cell2mat(paramplotdata(i,1:6)),'%g'),paramplotdata(i,7));
end

writetable(cell2table(paramplotdata,'VariableNames',{'model_data','pois_data','norm_data','model_fit','epidemic','informed','parameter','actual','estimate','uqid'}),...
    strcat(filename,'_paramplots','.csv'));
writetable(cell2table(fitplotdata,'VariableNames',{'generating_model','fitting_model','noise','starting','t','data'}),...
    strcat(filename,'_fitplots','.csv'));



%%%% Write data for LaTeX

latextable = {};
minitable = {};
andsizer = cell(length(latexparamnames)+2,1);
andsizer(:) = {'&'};
endsizer = cell(length(latexparamnames)+2,1);
endsizer(:) = {'\\ \hline'};

header = {'Parameters/AIC','&','&',modelnames{1},'&',modelnames{2},'&',modelnames{3},'&',modelnames{4},'&',modelnames{5},'\\ \hline'};

for i=1:length(datasets)
    minitable = [[latexparamnames; '$\Delta AIC$'; '$\Ro$'] andsizer];
    for j = 1:length(models)
        if j == 3
            miniparams = [paramests{i,j}(1:2) paramests{i,j}(5) paramests{i,j}(3)*0.2+paramests{i,j}(4)*0.8 paramests{i,j}(end)/0.2];
        else
            miniparams = [paramests{i,j}(1:4) paramests{i,j}(end)];
        end  
        minitable = [minitable andsizer [num2cell(round(miniparams(1:4),4)'); num2cell(round(miniparams(5),9)); aicfcn(datasets{i},paramests{i,j},gofs(i,j)); round(r0fcns{j}(paramests{i,j}),4)]];
    end
    minitable = [minitable endsizer];
    
    %Replace AICs with delta-AICs
    minAIC = min(cell2mat(minitable(end-1,[4,6,8,10,12])));
    for k=[4,6,8,10,12]
        minitable{end-1,k} = round(minitable{end-1,k} - minAIC);
    end
    
    latextable = [latextable; header; minitable];
end

header = {'Parameters/AIC','&','&',modelnames{1},'&',modelnames{2},'&',modelnames{3},'&',modelnames{4},'&',modelnames{5},'\\ \hline'};
latextable = [header; latextable];

%swap order so 2nd model is DR and 4th model is gamma, for both fitting and generating models
latextable(:,[6,10]) = latextable(:,[10,6]);
% latextable([9:15,23:29],:) = latextable([23:29,9:15],:);  %Before I put headers on each table
latextable([9:16,25:32],:) = latextable([25:32,9:16],:);    %With headers

writetable(cell2table(latextable),strcat(filename,'_latextable','.txt'),'Delimiter',' ');



%% Additional Calcs - FIMs for all models under baseline conditions
% % assuming no noise (i.e. FIM has no weighting matrix)
% 
% baseICdata = 0.01*50000;
% allranks = [];
% allCVs = {}; %cell array since they'll be different lengths with different numbers of parameters
% 
% for i=1:length(models)
%     FIM = fisher(tspan,x0fcns{i}(baseICdata,baseparams{i}(end)),baseparams{i},models{i},yfcn);
%     %CRCOV = inv(FIM);
%     allranks = [allranks rank(FIM,0.1)];
%     allCVs = [allCVs (sqrt(diag(inv(FIM)))./baseparams{i}')*100];
% end


