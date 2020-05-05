% rmpath(genpath(pwd))
addpath(genpath(pwd))
addpath(genpath('D:\Gmsh'))
addpath(genpath('F:\13.Code\Git\MatlabFunctions'))
addpath('D:\Mosek\8\toolbox\r2014a');

% test installation of Mosek
mosekdiag
% option = 'static';
option = 'kinematic';
prjRoot = 'F:\13.Code\Git\PlateAnalysis';
switch option
    case 'static'
        path_rm = [pwd '\5.LimitAnalysis\kinematicApproach'];
        rmpath(genpath(path_rm))
    case 'kinematic' 
        path_rm = [pwd '\5.LimitAnalysis\staticApproach'];
        rmpath(genpath(path_rm))
end
clc