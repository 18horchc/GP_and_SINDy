% Authors: Sara Amato, Andrea Arnold <anarnold@wpi.edu>
% Last updated: August 29, 2025
%
% Corresponds with the manuscript:
% 
% S. Amato and A. Arnold (2025) Data-driven modeling and prediction of 
% microglial cell dynamics in the ischemic penumbra. 
% Preprint available on arXiv: https://arxiv.org/abs/2404.10915


clear; close all; clc

%% Data
newtime1 = [0, 1, 2, 3, 7, 14]; %These are the points where we observe experimental data
datapointsM2 = [5, 78.33, 179.5, 126.4, 319, 136.67]'; %Experimental measurements of M2
datapointsM1 = [5, 27.5, 122.5, 139.8,  445, 816.67]'; %Experimental measurements of M1

t = [0:0.1:14]; %time that we want to fit the polynomial over

pM1four = polyfit(newtime1, datapointsM1, 4); %obtain coefficients for a fourth order polynomial fit to M1 data
polyM1four = polyval(pM1four, t); %Compute y values of polynomial for t

pM2quad = polyfit(newtime1, datapointsM2, 2); %obtain coefficients for a second order polynomial fit to M2 data
polyM2quad = polyval(pM2quad, t); %Compute y values of polynomial for t

newtime2 = [0:1:14]; %These are the time points that we will extract proxy data for
for i = 1:length(newtime2)
datapointsM21(i) = polyM2quad(1 + 10*(i-1)); %noiseless proxy data for M1
datapointsM11(i) = polyM1four(1 + 10*(i-1)); %noiseless proxy data for M1
end

noiselevel = 0.1; %set noise level
noise1 = normrnd(0, noiselevel*(std(datapointsM11)), size(datapointsM11));
noise2 = normrnd(0, noiselevel*(std(datapointsM21)), size(datapointsM21));

datapointsM12 = datapointsM11 + noise1; %Add noise to M1 proxy data

%Keep experimental measurements as they are
datapointsM12(1) = 5;
datapointsM12(2) = 27.5;
datapointsM12(3) = 122.5;
datapointsM12(4) = 139.8;
datapointsM12(8) = 445;
datapointsM12(15) = 816.67;

datapointsM22 = datapointsM21 + noise2; %Add noise to M2 proxy data

%Keep experimental measurements as they are
datapointsM22(1) = 5;
datapointsM22(2) = 78.33;
datapointsM22(3) = 179.5;
datapointsM22(4) = 126.4;
datapointsM22(8) = 319;
datapointsM22(15) = 136.67;

%% Plotting data

figure(1)
m = scatter(newtime2, datapointsM12, 'k');
m.Marker = 'o';
m.SizeData = 150;
hold on
s = scatter(newtime1, datapointsM1, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(t, polyM1four, 'k', 'LineWidth', 2.0)
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('M1(t)', 'fontsize',20)
set(gca, 'fontsize', 30)

figure(2)
m = scatter(newtime2, datapointsM22, 'r');
m.Marker = 'o';
m.SizeData = 150;
hold on
s = scatter(newtime1, datapointsM2, 'r', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(t, polyM2quad, 'r', 'LineWidth', 2.0)
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('M2(t)', 'fontsize',20)
set(gca, 'fontsize', 30)


%% ESINDy
N = 1000; %Number of ensemble members
newtime = t;
M1data = datapointsM12';
M2data = datapointsM22';

M1 = zeros(length(newtime), 1, N); %Set up M1 to store values for every ensemble member
M2 = zeros(length(newtime), 1, N); %Set up M2 to store values for every ensemble member

%Candidate functions
Theta = [ones(length(M1data)-1, 1) M1data(1:end-1) M2data(1:end-1) M1data(1:end-1).^2 M2data(1:end-1).^2 M1data(1:end-1).*M2data(1:end-1)];

X2 = [M1data(2:end) M2data(2:end)];
epsguessstored = zeros(width(Theta), 2, N);

%there are 6 different functions in our library
choices = [1:1:6];

estimatestored = zeros(length(newtime), 2, N);
for i = 1:N
    l = 6 - 1; %In each iteration, we use 5 candidate functions
    total_cols = 6; %There are 6 candidate functions in total
    p = randsample(total_cols, l); %pick 5 out of the 6
    p = sort(p); %sort them
   g = ismember(choices, p); %logical expression
            for j = 1:length(g) %this will be useful later for storing our coefficient guesses
            if g(j) == 1
            else
                stored(i) = j; 
            end
            end
clear g
for g = 1:l
    Thetait(:, g) = Theta(:, p(g)); %Create new candidate function matrix with only the randomly selected columns
    X2it = X2;
end
lambda = 0.01; %threshold value

%STLS
epsguess = Thetait\X2it; % initial guess: Least-squares

for c=1:10
smallinds = (abs(epsguess)<lambda); % find small coefficients
epsguess(smallinds)=0; % and threshold
biginds1 = ~smallinds(:,1); 
biginds2 = ~smallinds(:, 2);
% Regress dynamics onto remaining terms to find sparse Xi
epsguess(biginds1,1) = Thetait(:,biginds1)\X2it(:,1);
epsguess(biginds2,2) = Thetait(:,biginds2)\X2it(:,2);
end

%Store values for each iteration
epsguessstored(stored, :, i) = NaN; %If candidate function was not included, set to NaN

epsguessstored(1:stored(i)-1, :, i) = epsguess(1:stored(i) - 1, :); %store coeff values
epsguessstored(stored(i)+1:end, :, i) = epsguess(stored(i):end, :);

clear epsguess

end


%% Post-processing to get most likely model
%This is not an automated process and the user will have to sort through
%the models manually

epsguessstored(isnan(epsguessstored)) = 0; %Set NaN values to 0 for simplicity

indices = 1; %these are the ensemble numbers that we have identified
counter(1) = 1; %We count up how many iterations yield the same model
for i = 1:N
    if i ~= 1
    if epsguessstored(:, :, i) == epsguessstored(:, :, 1) 
        indices(end + 1) = i;
        counter(1) = counter(1) + 1;
    end
    end
end

%keep repeating this process until you've accounted for all models
counter(2) = 0;
for i = 1:N
    if epsguessstored(:, :, i) == epsguessstored(:, :, 2)
        indices(end + 1) = i;
        counter(2) = counter(2) + 1;
    end
end


%% Model Results
epsguess = epsguessstored(:, :, 1) %%%%%%Note: The user will have to change this index to be one of the most likely models

%Equations for M1 and M2
estimate = zeros(length(newtime2), 2);
estimate(1,:)=5;
for i = 2:(length(newtime2))
    estimate(i, 1) = epsguess(1, 1) + epsguess(2,1)*estimate(i-1, 1) + epsguess(3, 1)*estimate(i-1,2) + epsguess(4, 1)*estimate(i-1, 1)^2 + epsguess(5,1)*estimate(i-1, 2)^2 + epsguess(6, 1)*estimate(i-1, 1)*estimate(i-1, 2); 
    estimate(i, 2) = epsguess(1, 2) + epsguess(2,2)*estimate(i-1,1) + epsguess(3, 2)*estimate(i-1,2) + epsguess(4, 2)*estimate(i-1, 1)^2 + epsguess(5,2)*estimate(i-1, 2)^2 + epsguess(6, 2)*estimate(i-1, 1)*estimate(i-1, 2); 
end

%Plot
figure(201)
plot(newtime2, estimate(:, 1), 'k', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM1, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(newtime2, estimate(:, 2), 'r', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM2, 'r', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 30)

