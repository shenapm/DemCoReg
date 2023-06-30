clear; format shortG; warning off; close all; clc;
addpath('.\include\');

%% Project name
sPrj = 'HMA-2';

sFoldData = ['..\','Data\',sPrj,'\'];
sFoldRes  = ['..\','Results\',sPrj,'\']; mkdir(sFoldRes);

%% Method
cMethod = {...
    'N23: Nuth and Kääb standard version';...
    'N13: Nuth and Kääb simplified version';...
    'L23: Nuth and Kääb linear version';...
    'L57: Rosenholm and Torlegard (The revised version using terrain slope and aspect)';...
    'L57: Rosenholm and Torlegard (The raw version using terrain gradients)';};

%% Configuration parameters
% Colormap
colDiv = brewermap(256,'RdYlBu'); %Diverging 
colSeq = brewermap(256,'Blues');  %Sequential

% Regression
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
Nuth3 = @(b,x)(b(1)*cosd(b(2)-x)+b(3));  %Equation #3 in Nuth and Kääb (2011)
Nuth2 = @(b,x)(b(1)*cosd(b(2)-x(:,1)).*x(:,2)+b(3)); %Equation #2 in Nuth and Kääb (2011)

% Thresholds
tIterMaxNum = 20;       % Maximum number of iterations
tIterMadDif = 0.001;    % MAD difference
tMinAngS = 3.0;         % Minimum slope angle