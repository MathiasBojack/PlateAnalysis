clear
pathcontrol

%%
a  = 12; %m
b  = 120;
h  = 0.15;
%% plot the maximal deflection
TimeList = [60 1200:1200:14400];

figure(1)
clf;
pltGlobalSettings
figData.size = [600 200];
figData.xData = cell(1,1);
figData.yData = cell(1,1);
figData.plotType = 'normal';

for ii = 1: length(TimeList)
    Time = TimeList(ii);
    % KL
    fnameKL = [ prjRoot '\3.Thermoelastic\Kirchhoff-Love\DeformedShape\'...
        'Deformee_a_', num2str(a), '_b_', num2str(b) '_h_' ...
        num2str(h) '_Time_' num2str(Time) '_KL.mat'];
    KL = load(fnameKL);
    minKL(ii) = min(KL.Solution.W(:));
    % VK
    fnameVK = [ prjRoot '\3.Thermoelastic\vonKarman\DeformedShape\'...
        'Deformee_a_', num2str(a), '_b_', num2str(b) '_h_' ...
        num2str(h) '_Time_' num2str(Time) '.mat'];
    VK = load(fnameVK);
    minVK(ii) = min(VK.Solution.W(:));
    % lambda
%     delta(ii) = norm( (VK.Solution.W - KL.Solution.W))/ norm(KL.Solution.W);
    deltaMax(ii) = min(VK.Solution.W(:))/ min(KL.Solution.W(:)) -1;
end


% plot
pltIndex = 2:7;
figData.xData{1} = TimeList(pltIndex)'/3600;
figData.yData{1} = [minKL(pltIndex)' minVK(pltIndex)'];
figData.xlabel{1} = 'Time [h]';
figData.ylabel{1} = 'Deflection [m]';

figData.xData{2} = TimeList(pltIndex)'/3600;
figData.yData{2} = [ deltaMax(pltIndex)'];
figData.ylabel{2} = '$\delta$';


% legend
figData.legText{1} ='VK';
figData.legText{2} ='KL';
figData.legText{3} ='$\delta$';
% figData.Position =

[fig,figData] = multiAxisPlot(figData); % raw plot on different axes
[fig, figData] = reArrangeAxis(figData,fig);

%% figure (2) stability factor


figure(2)
clf;
pltGlobalSettings
figData.Number = 2;
figData.size = [400 400];
figData.xData = cell(1,1);
figData.yData = cell(1,1);
figData.plotType = 'normal';

% KL
fnameKL = [ prjRoot '\7.CaseStudy\4sPlate\KLmodel\'...
    'lambda_a_', num2str(a), '_b_', num2str(b) '_h_' ...
    num2str(h) '_KL.mat'];
KL = load(fnameKL);

% VK
fnameVK = [ prjRoot '\7.CaseStudy\4sPlate\VKmodel\'...
    'lambda_a_', num2str(a), '_b_', num2str(b) '_h_' ...
    num2str(h) '_VK.mat'];
VK = load(fnameVK);

% plot

figData.xData{1} = KL.Time_list(pltIndex)'/3600;
figData.yData{1} = [KL.lamb_stat(pltIndex) VK.lamb_stat(pltIndex)];
figData.xlabel{1} = 'Time [h]';
figData.ylabel{1} = '$\lambda$';
figData.legText{1} ='VK';
figData.legText{2} ='KL';

[fig,figData] = multiAxisPlot(figData); % raw plot on different axes
[fig, figData] = reArrangeAxis(figData,fig);