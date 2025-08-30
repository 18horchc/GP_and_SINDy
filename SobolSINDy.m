% Authors: Sara Amato, Andrea Arnold <anarnold@wpi.edu>
% Last updated: August 29, 2025
%
% Corresponds with the manuscript:
% 
% S. Amato and A. Arnold (2025) Data-driven modeling and prediction of 
% microglial cell dynamics in the ischemic penumbra. 
% Preprint available on arXiv: https://arxiv.org/abs/2404.10915


clear; close all; clc;

tic %set timer

%% Sobol method

%Create M by p sample matrices A and B
M=25000; %Set M
p= 6; %P is the number of parameters

% Generate sobol sequences
ss = sobolset(2*p, 'Skip', 1);
X = net(ss, M);

q1 = X(:, 1);
q1hat = X(:, p+1);
q2 = X(:, 2);
q2hat = X(:, p+2);
q3 = X(:, 3);
q3hat = X(:, p+3);
q4 = X(:, 4);
q4hat = X(:, p+4);
q5 = X(:, 5);
q5hat = X(:, p+5);
q6 = X(:, 6);
q6hat = X(:, p+6);

%SINDy point estimates
nonzerocoeffM1 = [32.9650, -0.1377, 0.3157];
nonzerocoeffM2 = [70.3415, -0.1569, 0.0217];

%Rescale 0.8 to 1.2 point estimate value
%Parameter: (1)
p1_low = 0.80*nonzerocoeffM1(1); 
p1_high = 1.2*nonzerocoeffM1(1);
q1 = p1_low + (p1_high-p1_low).*q1;
q1hat = p1_low + (p1_high-p1_low).*q1hat;

%Parameter: (2)
p2_low = 0.80*nonzerocoeffM1(2); 
p2_high = 1.2*nonzerocoeffM1(2);
q2 = p2_low + (p2_high - p2_low).*q2;
q2hat = p2_low + (p2_high - p2_low).*q2hat;

%Parameter: (3)
p3_low =1.2*nonzerocoeffM1(3); 
p3_high = 0.8*nonzerocoeffM1(3);
q3 = p3_low + (p3_high - p3_low).*q3;
q3hat = p3_low + (p3_high - p3_low).*q3hat;

%Parameter: (4)
p4_low = 0.80*nonzerocoeffM2(1); 
p4_high = 1.2*nonzerocoeffM2(1);
q4 = p4_low + (p4_high - p4_low).*q4;
q4hat = p4_low + (p4_high - p4_low).*q4hat;

%Parameter: (5)
p5_low = 0.80*nonzerocoeffM2(2); 
p5_high = 1.2*nonzerocoeffM2(2);
q5 = p5_low + (p5_high - p5_low).*q5;
q5hat = p5_low + (p5_high - p5_low).*q5hat;

%Parameter:  (6)
p6_low = 0.80*nonzerocoeffM2(3); 
p6_high = 1.2*nonzerocoeffM2(3);
q6 = p6_low + (p6_high - p6_low).*q6;
q6hat = p6_low + (p6_high - p6_low).*q6hat;


%create matrices for Sobol Sensitivity 
A = [q1, q2, q3, q4, q5, q6]; 
B = [q1hat, q2hat, q3hat, q4hat, q5hat, q6hat];

C1 = [q1, q2hat, q3hat, q4hat, q5hat, q6hat];
C2 = [q1hat, q2, q3hat, q4hat, q5hat, q6hat];
C3 = [q1hat, q2hat, q3, q4hat, q5hat, q6hat];
C4 = [q1hat, q2hat, q3hat, q4, q5hat, q6hat];
C5 = [q1hat, q2hat, q3hat, q4hat, q5, q6hat];
C6 = [q1hat, q2hat, q3hat, q4hat, q5hat, q6]; 


%Compute M by 1 outputs of model vectors
for j=1:M
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], A(j, :));
    yAM1(j) = trapint(time, sol(:, 1)); 
    yAM2(j) = trapint(time, sol(:, 2)); 
    yAratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yAsum(j) = trapint(time, sol(:, 1)+sol(:, 2)); 
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], B(j, :));
    yBM1(j) = trapint(time, sol(:, 1));   
    yBM2(j) = trapint(time, sol(:, 2)); 
    yBratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yBsum(j) = trapint(time, sol(:, 1)+sol(:, 2)); 
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C1(j, :));
    yC1M1(j) = trapint(time, sol(:, 1)); 
    yC1M2(j) = trapint(time, sol(:, 2)); 
    yC1ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC1sum(j) = trapint(time, sol(:, 1)+sol(:, 2)); 
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C2(j, :));
    yC2M1(j) = trapint(time, sol(:, 1)); 
    yC2M2(j) = trapint(time, sol(:, 2)); 
    yC2ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC2sum(j) = trapint(time, sol(:, 1)+sol(:, 2));
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C3(j, :));
    yC3M1(j) = trapint(time,  sol(:, 1)); 
    yC3M2(j) = trapint(time, sol(:, 2)); 
    yC3ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC3sum(j) = trapint(time, sol(:, 1)+sol(:, 2));
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C4(j, :));
    yC4M1(j) = trapint(time,  sol(:, 1));
    yC4M2(j) = trapint(time, sol(:, 2)); 
    yC4ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC4sum(j) = trapint(time, sol(:, 1)+sol(:, 2));
    
    clear sol
    
    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C5(j, :));
    yC5M1(j) = trapint(time, sol(:, 1));
    yC5M2(j) = trapint(time, sol(:, 2)); 
    yC5ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC5sum(j) = trapint(time, sol(:, 1)+sol(:, 2));
    
    clear sol

    [time, sol] = ode15s(@rhs_DESINDy, 0:0.1:14, [5,5], [], C6(j, :));
    yC6M1(j) = trapint(time, sol(:, 1));
    yC6M2(j) = trapint(time, sol(:, 2)); 
    yC6ratio(j) = trapint(time, sol(:, 1)./sol(:, 2)); 
    yC6sum(j) = trapint(time, sol(:, 1)+sol(:, 2));
  
end

%% Estimate first order sesnsitivities
f02M1 = ((1/M)*sum(yAM1))*((1/M)*sum(yBM1));
denominatorM1 = ((1/M)*(yAM1*yAM1')- (f02M1));

f02M2 = ((1/M)*sum(yAM2))*((1/M)*sum(yBM2));
denominatorM2 = ((1/M)*(yAM2*yAM2')- (f02M2));

f02ratio = ((1/M)*sum(yAratio))*((1/M)*sum(yBratio));
denominatorratio = ((1/M)*(yAratio*yAratio')- (f02ratio));

f02sum = ((1/M)*sum(yAsum))*((1/M)*sum(yBsum));
denominatorsum = ((1/M)*(yAsum*yAsum')- (f02sum));

%M1
S1M1 = ((((1/M)*(yAM1*yC1M1'))-(f02M1))/denominatorM1);
S2M1 = ((((1/M)*(yAM1*yC2M1'))-(f02M1))/denominatorM1);
S3M1 = ((((1/M)*(yAM1*yC3M1'))-(f02M1))/denominatorM1);
S4M1 = ((((1/M)*(yAM1*yC4M1'))-(f02M1))/denominatorM1);
S5M1 = ((((1/M)*(yAM1*yC5M1'))-(f02M1))/denominatorM1);
S6M1 = ((((1/M)*(yAM1*yC6M1'))-(f02M1))/denominatorM1);


sens_vec = [S1M1, S2M1, S3M1, S4M1, S5M1, S6M1]; 
[rank_Sob rank_Sobp] = sort(sens_vec, 'descend');


%M2
S1M2 = ((((1/M)*(yAM2*yC1M2'))-(f02M2))/denominatorM2);
S2M2 = ((((1/M)*(yAM2*yC2M2'))-(f02M2))/denominatorM2);
S3M2 = ((((1/M)*(yAM2*yC3M2'))-(f02M2))/denominatorM2);
S4M2 = ((((1/M)*(yAM2*yC4M2'))-(f02M2))/denominatorM2);
S5M2 = ((((1/M)*(yAM2*yC5M2'))-(f02M2))/denominatorM2);
S6M2 = ((((1/M)*(yAM2*yC6M2'))-(f02M2))/denominatorM2);


sens_vecM2 = [S1M2, S2M2, S3M2, S4M2, S5M2, S6M2];
[rank_SobM2 rank_SobpM2] = sort(sens_vecM2, 'descend');

%% Compute Total Sensitivities

%M1
ST1M1 =  1- ((((1/M)*(yBM1*yC1M1')) - (f02M1))/denominatorM1);
ST2M1 =  1- ((((1/M)*(yBM1*yC2M1')) - (f02M1))/denominatorM1);
ST3M1 =  1- ((((1/M)*(yBM1*yC3M1')) - (f02M1))/denominatorM1);
ST4M1 =  1- ((((1/M)*(yBM1*yC4M1')) - (f02M1))/denominatorM1);
ST5M1 =  1- ((((1/M)*(yBM1*yC5M1')) - (f02M1))/denominatorM1);
ST6M1 =  1- ((((1/M)*(yBM1*yC6M1')) - (f02M1))/denominatorM1);

STM1 = [ST1M1 ST2M1 ST3M1, ST4M1, ST5M1, ST6M1];

%M2
ST1M2 =  1- ((((1/M)*(yBM2*yC1M2')) - (f02M2))/denominatorM2);
ST2M2 =  1- ((((1/M)*(yBM2*yC2M2')) - (f02M2))/denominatorM2);
ST3M2 =  1- ((((1/M)*(yBM2*yC3M2')) - (f02M2))/denominatorM2);
ST4M2 =  1- ((((1/M)*(yBM2*yC4M2')) - (f02M2))/denominatorM2);
ST5M2 =  1- ((((1/M)*(yBM2*yC5M2')) - (f02M2))/denominatorM2);
ST6M2 =  1- ((((1/M)*(yBM2*yC6M2')) - (f02M2))/denominatorM2);

STM2 = [ST1M2 ST2M2 ST3M2, ST4M2, ST5M2, ST6M2];


%% Generalized (this is what is used in the present work)

STgen1 = (ST1M1*denominatorM1 + ST1M2*denominatorM2)/(denominatorM1 + denominatorM2);
STgen2 = (ST2M1*denominatorM1 + ST2M2*denominatorM2)/(denominatorM1 + denominatorM2);
STgen3 = (ST3M1*denominatorM1 + ST3M2*denominatorM2)/(denominatorM1 + denominatorM2);
STgen4 = (ST4M1*denominatorM1 + ST4M2*denominatorM2)/(denominatorM1 + denominatorM2);
STgen5 = (ST5M1*denominatorM1 + ST5M2*denominatorM2)/(denominatorM1 + denominatorM2);
ST6gen = (ST6M1*denominatorM1 + ST6M2*denominatorM2)/(denominatorM1 + denominatorM2);

STgen = [STgen1, STgen2, STgen3, STgen4, STgen5, ST6gen];


toc %end timer

%% Plot- gives figure 7

figure(2)
labels = [1,2,3, 4, 5, 6]; 
plot([1:1:6], STgen,  'or', 'MarkerSize',15, 'MarkerEdgeColor','blue','MarkerFaceColor',[0 0 1])
hold on
xticks(labels)
ylim([0 0.8])
xlabel('Parameter Indices')
ylabel('$\tilde{S}_{T_i}$', 'Interpreter','latex')
set(gca, 'fontsize', 24)


%% Function

function rhs = rhs_DESINDy(t, inits, params)
M1 = inits(1);
M2 = inits(2);

theta1 = params(1);
theta2 = params(2);
theta3 = params(3);
theta4 = 0;
theta5 = 0;
theta6 = 0;
theta7 = params(4);
theta8 = params(5);
theta9 = params(6);
theta10 = 0;
theta11 = 0;
theta12 = 0;

% Equations
dM1 = theta1 + theta2*M1 + theta3*M2 + theta4*M1^2 + theta5*M2^2 + theta6*M1*M2;
dM2 = theta7 + theta8*M1 + theta9*M2 + theta10*M1^2 + theta11*M2^2 + theta12*M1*M2;

rhs = [dM1; dM2];
end

