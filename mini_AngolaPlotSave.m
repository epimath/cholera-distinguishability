% This file plots and saves data for fitting to Angola
% 1/14/2016 Marisa Eisenberg (marisae@umich.edu)

% Run this after the fitting section in SimFitAllMain.m

% This is basically a slightly modified version of the plot & save sections
% in SimFitAllMain.m. At some point I should really make this nicer so it's
% all one thing with options depending on what you're doing. 

%% Plot all fits

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
    title(strcat('Angola Data'))
    ylabel('Cases');
    xlabel('Days');
    legend(['Data' modelnames]);
    
    %optional save figs
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    saveas(fig,strcat(filename,'_','Angola Data'),'pdf')
    saveas(fig,strcat(filename,'_','Angola Data'),'fig')
end

%% AICs

AICs = []
for j = 1:length(models)
AICs = [AICs aicfcn(datasets{1},paramests{1,j},gofs(1,j))];
end


%% Save data


%%%% Matlab Data
save(strcat(filename,'.mat'))



%% Write data for LaTeX

latextable = {};
minitable = {};
andsizer = cell(length(latexparamnames)+2,1);
andsizer(:) = {'&'};
endsizer = cell(length(latexparamnames)+2,1);
endsizer(:) = {'\\ \hline'};

header = {'Parameters/AIC','&','&',modelnames{1},'&',modelnames{2},'&',modelnames{3},'&',modelnames{4},'&',modelnames{5},'\\ \hline'};

for i=1:length(datasets)
%     minitable = [[latexparamnames; '$\Delta AIC$'; '$\Ro$'] andsizer];
    minitable = [[latexparamnames; '$AIC$'; '$\Ro$'] andsizer];
    for j = 1:length(models)
        if j == 3
            miniparams = [paramests{i,j}(1:2) paramests{i,j}(5) paramests{i,j}(3)*0.2+paramests{i,j}(4)*0.8 paramests{i,j}(end)/0.2];
        else
            miniparams = [paramests{i,j}(1:4) paramests{i,j}(end)];
        end  
        minitable = [minitable andsizer [num2cell(round(miniparams(1:4),4)'); num2cell(round(miniparams(5),9)); round(aicfcn(datasets{i},paramests{i,j},gofs(i,j))); round(r0fcns{j}(paramests{i,j}),4)]];
    end
    minitable = [minitable endsizer];
    
    %Replace AICs with delta-AICs
%     minAIC = min(cell2mat(minitable(end-1,[4,6,8,10,12])));
%     for k=[4,6,8,10,12]
%         minitable{end-1,k} = round(minitable{end-1,k} - minAIC);
%     end
    
    latextable = [latextable; header; minitable];
end

header = {'Parameters/AIC','&','&',modelnames{1},'&',modelnames{2},'&',modelnames{3},'&',modelnames{4},'&',modelnames{5},'\\ \hline'};
latextable = [header; latextable];

%swap order so 2nd model is DR and 4th model is gamma, for both fitting and generating models
latextable(:,[6,10]) = latextable(:,[10,6]);
% latextable([9:15,23:29],:) = latextable([23:29,9:15],:);  %Before I put headers on each table
% latextable([9:16,25:32],:) = latextable([25:32,9:16],:);    %With headers

writetable(cell2table(latextable),strcat(filename,'_latextable','.txt'),'Delimiter',' ');



%% Write data for fit plots
paramplotdata = {};
fitplotdata = {};

psizer = ones(length(dataparams{1}),1);
dsizer = cell(size(tspan));
msizer = cell(size(tspan));
nssizer = {};
for i=1:length(tspan)
    nssizer = [nssizer; noisestart];
end

% for i=1:length(datasets)
    for j = 1:length(models)
%         if j == 3
%             miniparams = [paramests{i,j}(1:2) paramests{i,j}(5) paramests{i,j}(3)*0.2+paramests{i,j}(4)*0.8 paramests{i,j}(end)/0.2];
%         else
%             miniparams = [paramests{i,j}(1:4) paramests{i,j}(end)];
%         end  
%         paramplotdata = [paramplotdata; num2cell([i*psizer extracols(1)*psizer extracols(2)*psizer j*psizer extracols(3)*psizer extracols(4)*psizer]) paramnames num2cell([dataparams{1}' miniparams'])];
        
        dsizer(:) = {'Angola'}; %modelnames(i);
        msizer(:) = modelnames(j);
        fitplotdata = [fitplotdata; dsizer msizer nssizer num2cell(tspan) num2cell(fits{1,j})];
    end
    
    % add true data to fitplot data
    dsizer(:) = {'Angola'}; %modelnames(i);
    msizer(:) = {'True Data'};
    fitplotdata = [fitplotdata; dsizer msizer nssizer num2cell(tspan) num2cell(datasets{1})];
% end


%    write the little labels at the end - could merge this with above, but meh
% for i=1:size(paramplotdata,1)
%     paramplotdata{i,10} = strcat(num2str(cell2mat(paramplotdata(i,1:6)),'%g'),paramplotdata(i,7));
% end

% writetable(cell2table(paramplotdata,'VariableNames',{'model_data','pois_data','norm_data','model_fit','epidemic','informed','parameter','actual','estimate','uqid'}),...
%     strcat(filename,'_paramplots','.csv'));
writetable(cell2table(fitplotdata,'VariableNames',{'generating_model','fitting_model','noise','starting','t','data'}),...
    strcat(filename,'_fitplots','.csv'));
