% Authors: Sara Amato, Andrea Arnold <anarnold@wpi.edu>
% Last updated: August 29, 2025
%
% Corresponds with the manuscript:
% 
% S. Amato and A. Arnold (2025) Data-driven modeling and prediction of 
% microglial cell dynamics in the ischemic penumbra. 
% Preprint available on arXiv: https://arxiv.org/abs/2404.10915

clear; close all; clc

% Plots- for each section for ESINDy

%% Average data plot for all 10% noise levels- this will give a plot similar to figure 12
%resulting epsguess parameter gives parameters in (3.4) and (3.5)

load  DBN10noisedatatry1
epsguess1 = epsguess;

clear epsguess

load  DBN10noisedatatry2
epsguess2 = epsguess;

clear epsguess

load  DBN10noisedatatry3
epsguess3 = epsguess;

clear epsguess

load DBN10noisedatatry4
epsguess4 = epsguess;

clear epsguess

load  DBN10noisedatatry5
epsguess5 = epsguess;

clear epsguess

load  DBN10noisedatatry6
epsguess6 = epsguess;

clear epsguess

load  DBN10noisedatatry7
epsguess7 = epsguess;

clear epsguess

load  DBN10noisedatatry8
epsguess8 = epsguess;

clear epsguess

load  DBN10noisedatatry9
epsguess9 = epsguess;

clear epsguess

load  DBN10noisedatatry10
epsguess10 = epsguess;


A(:, :, 1) = epsguess1;
A(:, :, 2) = epsguess2;
A(:, :, 3) = epsguess3;
A(:, :, 4) = epsguess4;
A(:, :, 5) = epsguess5;
A(:, :, 6) = epsguess6;
A(:, :, 7) = epsguess7;
A(:, :, 8) = epsguess8;
A(:, :, 9) = epsguess9;
A(:, :, 10) = epsguess10;


epsguess = mean(A, 3)

newtime2 = [0:1:14];
estimate = zeros(length(newtime2), 2);
estimate(1,:)=5;
for i = 2:(length(newtime2))
    estimate(i, 1) = epsguess(1, 1) + epsguess(2,1)*estimate(i-1, 1) + epsguess(3, 1)*estimate(i-1,2) + epsguess(4, 1)*estimate(i-1, 1)^2 + epsguess(5,1)*estimate(i-1, 2)^2 + epsguess(6, 1)*estimate(i-1, 1)*estimate(i-1, 2); %+epsguess(7,1)*estimate(i-1, 1)^3 + epsguess(8,1)*estimate(i-1, 2)^3;
    estimate(i, 2) = epsguess(1, 2) + epsguess(2,2)*estimate(i-1,1) + epsguess(3, 2)*estimate(i-1,2) + epsguess(4, 2)*estimate(i-1, 1)^2 + epsguess(5,2)*estimate(i-1, 2)^2 + epsguess(6, 2)*estimate(i-1, 1)*estimate(i-1, 2); %+epsguess(7,1)*estimate(i-1, 1)^3 + epsguess(8,1)*estimate(i-1, 2)^3;
end


newtime1 = [0, 1, 2, 3, 7, 14];
datapointsM2real = [5, 78.33, 179.5, 126.4, 319, 136.67]';
datapointsM1real = [5, 27.5, 122.5, 139.8,  445, 816.67]';
figure(201)
plot(newtime2, estimate(:, 1), 'k', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM1real, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(newtime2, estimate(:, 2), 'r', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM2real, 'r', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 30)

%% This is the plot for uncertainty forecasts and predictions
%Figure 13 in the paper

%Run this one after the average and don't clear all... before running this
%section
newtime3 = [0:1:50];
estimate_pred_pts_M1 = [estimate(1,1); estimate(2,1); estimate(3,1); estimate(4,1); estimate(8,1); estimate(15,1)];
RMSE_for_predM1 = zeros(1, length(newtime3));
RMSE_for_predM1(1) = rmse(estimate_pred_pts_M1(1), datapointsM1real(1));
RMSE_for_predM1(2) = rmse(estimate_pred_pts_M1(2), datapointsM1real(2));
RMSE_for_predM1(3) = rmse(estimate_pred_pts_M1(3), datapointsM1real(3));
RMSE_for_predM1(4:7) = rmse(estimate_pred_pts_M1(4), datapointsM1real(4));
RMSE_for_predM1(8:14) = rmse(estimate_pred_pts_M1(5), datapointsM1real(5));
RMSE_for_predM1(15:end) = rmse(estimate_pred_pts_M1(6), datapointsM1real(6));

estimate_pred_pts_M2 = [estimate(1,2); estimate(2,2); estimate(3,2); estimate(4,2); estimate(8,2); estimate(15,2)];
RMSE_for_predM2 = zeros(1, length(newtime3));
RMSE_for_predM2(1) = rmse(estimate_pred_pts_M2(1), datapointsM2real(1));
RMSE_for_predM2(2) = rmse(estimate_pred_pts_M2(2), datapointsM2real(2));
RMSE_for_predM2(3) = rmse(estimate_pred_pts_M2(3), datapointsM2real(3));
RMSE_for_predM2(4:7) = rmse(estimate_pred_pts_M2(4), datapointsM2real(4));
RMSE_for_predM2(8:14) = rmse(estimate_pred_pts_M2(5), datapointsM2real(5));
RMSE_for_predM2(15:end) = rmse(estimate_pred_pts_M2(6), datapointsM2real(6));

estimate1 = zeros(length(newtime3), 2);
estimate1(1,:)=5;
for i = 2:(length(newtime3))
    estimate1(i, 1) = epsguess(1, 1) + epsguess(2,1)*estimate1(i-1, 1) + epsguess(3, 1)*estimate1(i-1,2) + epsguess(4, 1)*estimate1(i-1, 1)^2 + epsguess(5,1)*estimate1(i-1, 2)^2 + epsguess(6, 1)*estimate1(i-1, 1)*estimate1(i-1, 2); %+epsguess(7,1)*estimate(i-1, 1)^3 + epsguess(8,1)*estimate(i-1, 2)^3;
    estimate1(i, 2) = epsguess(1, 2) + epsguess(2,2)*estimate1(i-1,1) + epsguess(3, 2)*estimate1(i-1,2) + epsguess(4, 2)*estimate1(i-1, 1)^2 + epsguess(5,2)*estimate1(i-1, 2)^2 + epsguess(6, 2)*estimate1(i-1, 1)*estimate1(i-1, 2); %+epsguess(7,1)*estimate(i-1, 1)^3 + epsguess(8,1)*estimate(i-1, 2)^3;
end

for j = 1:length(newtime3)
    CIM1mean(j) = estimate1(j, 1);
    CIM1std(j) = RMSE_for_predM1(j);
    z = 1.96;
    upperCIM1(j) = CIM1mean(j) + z*CIM1std(j);
    lowerCIM1(j) = CIM1mean(j) - z*CIM1std(j);
    
    CIM2mean(j) = estimate1(j, 2);
    CIM2std(j) = RMSE_for_predM2(j);
    z = 1.96;
    upperCIM2(j) = CIM2mean(j) + z*CIM2std(j);
    lowerCIM2(j) = CIM2mean(j) - z*CIM2std(j);
    
end

figure(202)
x = newtime3;
fill( [x fliplr(x)],  [upperCIM1 fliplr(lowerCIM1)], 'k','EdgeColor','none','FaceAlpha',0.2);
hold on
fill( [x fliplr(x)],  [upperCIM2 fliplr(lowerCIM2)],'r','EdgeColor','none', 'FaceAlpha',0.2);
hold on
plot(newtime3, estimate1(:, 1), 'k', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM1real, 'k', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
plot(newtime3, estimate1(:, 2), 'r', 'LineWidth',2.0)
hold on
s = scatter(newtime1, datapointsM2real, 'r', 'filled');
s.Marker = 'hexagram';
s.SizeData = 150;
hold on
scatter(35, 550,  'k', 'filled',  'SizeData', 150)
hold on
scatter(35, 50,  'r', 'filled',  'SizeData', 150)
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
xlim([0 50])
ylabel('cells/mm^2', 'fontsize',20)
set(gca, 'fontsize', 30)


%% Plot for proxy data - figure 3

clear; close all; clc

load DBNnoisytry2.mat

datapointsM1four = [polyM1four(41); polyM1four(51); polyM1four(61); ...
    polyM1four(81); polyM1four(91); polyM1four(101); ...
    polyM1four(111); polyM1four(121); polyM1four(131)]; 

datapointsM2quad = [polyM2quad(41); polyM2quad(51); polyM2quad(61); ...
    polyM2quad(81); polyM2quad(91); polyM2quad(101); ...
    polyM2quad(111); polyM2quad(121); polyM2quad(131)]; 

newtimefive = [4, 5, 6, 8, 9, 10, 11, 12, 13];
errM1 = 2*noiselevel*(std(datapointsM11))*ones(size(newtimefive));
errM2 = 2*noiselevel*(std(datapointsM21))*ones(size(newtimefive));

errM10 = 2*0.5*std([0, 10]);
errM11 = 2*0.5*std([5, 50]);
errM12 = 2*0.5*std([120, 125]);
errM13 = 2*0.5*std([375, 62, 102, 100, 60]);
errM17 = 2*0.5*std([600, 55, 900, 225]);
errM114 = 2*0.5*std([750, 400, 1300]);

figure(999)
s = scatter(newtime1, datapointsM1real, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
plot(t, polyM1four, 'k', 'LineWidth', 2.0)
hold on
errorbar(newtimefive,datapointsM1four,errM1, 'k', "LineStyle", "none")
hold on
s = scatter(0, 0, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(0, 10, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(1, 5, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(1, 50, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(2, 120, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(2, 125, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 375, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 62, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 102, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 100, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 60, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 600, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 55, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 900, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 225, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14, 750, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14, 400, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14, 1300, 'k', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 1500])
ylabel('M1(t)', 'fontsize',20)
set(gca, 'fontsize', 30)


figure(9999)
s = scatter(newtime1, datapointsM2real, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
plot(t, polyM2quad, 'r', 'LineWidth', 2.0)
hold on
errorbar(newtimefive,datapointsM2quad,errM2, 'r', 'LineStyle', "none")
hold on
s = scatter(0, 0, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(0, 10, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(1, 170, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(1, 15, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(1, 50, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(2, 90, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(2, 269, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 300, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 15, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 57, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 100, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(3, 160, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 600, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7,6, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 400, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(7, 270, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14,200, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14, 110, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold on
s = scatter(14, 100, 'r', 'filled');
s.Marker = 'o';
s.SizeData = 10;
hold off
xlabel('Time (Days)', 'fontsize',20) 
ylim([0 500])
ylabel('M2(t)', 'fontsize',20)
set(gca, 'fontsize', 30)


%% Plots for parameters- Figure 11

clear; close all; clc;

load  DBN10noisedatatry1
epsguess1 = epsguess;

clear epsguess

load  DBN10noisedatatry2
epsguess2 = epsguess;

clear epsguess

load  DBN10noisedatatry3
epsguess3 = epsguess;

clear epsguess

load  DBN10noisedatatry4
epsguess4 = epsguess;

clear epsguess

load  DBN10noisedatatry5
epsguess5 = epsguess;

clear epsguess

load  DBN10noisedatatry6
epsguess6 = epsguess;

clear epsguess

load  DBN10noisedatatry7
epsguess7 = epsguess;

clear epsguess

load  DBN10noisedatatry8
epsguess8 = epsguess;

clear epsguess

load  DBN10noisedatatry9
epsguess9 = epsguess;

clear epsguess

load  DBN10noisedatatry10
epsguess10 = epsguess;


A(:, :, 1) = epsguess1;
A(:, :, 2) = epsguess2;
A(:, :, 3) = epsguess3;
A(:, :, 4) = epsguess4;
A(:, :, 5) = epsguess5;
A(:, :, 6) = epsguess6;
A(:, :, 7) = epsguess7;
A(:, :, 8) = epsguess8;
A(:, :, 9) = epsguess9;
A(:, :, 10) = epsguess10;


coeff1 = A(1, 1, :);
coeff1 = reshape(coeff1, [], 1);
coeff2 = A(1, 2, :);
coeff2 = reshape(coeff2, [], 1);
coeff3 = A(2, 1, :);
coeff3 = reshape(coeff3, [], 1);
coeff4 = A(2, 2, :);
coeff4 = reshape(coeff4, [], 1);
coeff5 = A(3, 1, :);
coeff5 = reshape(coeff5, [], 1);
coeff6 = A(3, 2, :);
coeff6 = reshape(coeff6, [], 1);
figure(10)
boxplot([coeff1, coeff2],'Labels',{'$\theta_1$','$\theta_4$'})
set(gca,'TickLabelInterpreter','latex');
set(gca, 'FontSize', 30)


figure(11)
boxplot([coeff3, coeff5, coeff4, coeff6],'Labels',{'$\theta_2$','$\theta_3$','$\theta_5$', '$\theta_6$'})
set(gca,'TickLabelInterpreter','latex');
set(gca, 'FontSize', 30)

