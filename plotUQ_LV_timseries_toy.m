function plotUQ_LV_timeseries_toy(Xi,Xibs,Xistd,x0,tspan,polys,nUQ,pct,nE,tobs,xobs,options,lhpop)
% Rewritten for toy problem stability and dynamic scaling

rng(1,'twister')
iEE = 1;
che = 1000; % Increased from 10 to 1000 so toy populations (starting at 10) aren't ignored

% Pre-allocate for ensemble trajectories
numTime = length(tspan);
estimateUQ1 = zeros(numTime, nUQ);
estimateUQ2 = zeros(numTime, nUQ);

if nE == 1
    for iE = 1:nUQ
        XiUQ = Xibs(:,:,iE);
        [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);
        if max(max(xSINDYb)) < che
            estimateUQ1(:,iEE) = xSINDYb(:,1);
            estimateUQ2(:,iEE) = xSINDYb(:,2);
            iEE = iEE + 1;
        end
    end
else
    for iE = 1:nUQ
        iMod = randi(size(Xibs,3),1,nE);
        XiUQ = mean(Xibs(:,:,iMod),3);
        [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);
        if max(max(xSINDYb)) < che
            estimateUQ1(:,iEE) = xSINDYb(:,1);
            estimateUQ2(:,iEE) = xSINDYb(:,2);
            iEE = iEE + 1;
        end
    end
end

% Cleanup empty entries if some simulations diverged
estimateUQ1(:,iEE:end) = [];
estimateUQ2(:,iEE:end) = [];

% Ensemble Mean Simulation
[tspanSINDyXi,xSINDYXi]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polys),tspan,x0,options);

pctL = (100-pct)/2;
pctH = (100+pct)/2;
prcTT1 = prctile(estimateUQ1,[pctL 50 pctH],2);
prcTT2 = prctile(estimateUQ2,[pctL 50 pctH],2);

%% Plotting Logic
C10 = [0 128 255]/255; % Lynx Blue
CHare = [240,59,32]./255; % Hare Red
cGrey = 0.6*[1 1 1];
faceAlpha = 0.2; ms = 8; lw = 1.0; fos = 12;

% TOY PROBLEM FIX: Removed '+ 1900' time shift
tTrue = tobs; 
tspanPlot = tspan;

figure('Position', [100 100 700 450])

% Dynamic Scaling: Use actual std of toy data rows (3 and 2)
% This ensures 'x' markers align with the shaded clouds
n1 = std(lhpop(3,:)'); 
n2 = std(lhpop(2,:)');

% Hare Plotting
patch([tspanPlot'; flipud(tspanPlot')], [prcTT1(:,1)*n1; flipud(prcTT1(:,3)*n1)], CHare, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue, xobs(:,1)*n1, 'x', 'Color', CHare, 'LineWidth', lw, 'MarkerSize', ms); hold on
plot(tspanSINDyXi, xSINDYXi(:,1)*n1, '-', 'Color', CHare, 'LineWidth', 1.5); hold on

% Lynx Plotting
patch([tspanPlot'; flipud(tspanPlot')], [prcTT2(:,1)*n2; flipud(prcTT2(:,3)*n2)], C10, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue, xobs(:,2)*n2, 'x', 'Color', C10, 'LineWidth', lw, 'MarkerSize', ms); hold on
plot(tspanSINDyXi, xSINDYXi(:,2)*n2, '-', 'Color', C10, 'LineWidth', 1.5); hold on

% Formatting Fixes for 0-20 Range
ylabel('Population','interpreter','latex','FontSize',fos)
xlabel('Time','interpreter','latex','FontSize',fos)
box on; grid on;
legend({'Hare 95%','Hare Obs','Hare Model','Lynx 95%','Lynx Obs','Lynx Model'},'Location','northeastoutside');
title('Toy Problem: SINDy Uncertainty Quantification');