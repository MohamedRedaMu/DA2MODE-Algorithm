clc;
clear;
close all;

%% User settings
CECyear = 2020;   % Supported: 2017, 2020, 2022
fNo     = 3;      % Function number
nd      = 20;     % Dimension

lb = -100;        % Lower bound
ub = 100;         % Upper bound

%% Run DA2MODE
[goalReached, GlobalBest, countFE] = DA2MODE_algorithm(CECyear, fNo, nd, lb, ub);

%% Display results
disp('========================================');
disp('Algorithm            : DA2MODE');
disp(['CEC Year             : ', num2str(CECyear)]);
disp(['Function Number      : ', num2str(fNo)]);
disp(['Dimension            : ', num2str(nd)]);
disp(['Goal Reached         : ', num2str(goalReached)]);
disp(['Best Cost            : ', num2str(GlobalBest.Cost, '%.12g')]);
disp(['Function Evaluations : ', num2str(countFE)]);
disp('Best Position:');
disp(GlobalBest.Position);
disp('========================================');