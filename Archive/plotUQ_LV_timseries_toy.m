function plotUQ_LV_timseries_toy(Xi,Xibs,Xistd,x0,tspan,polys,nUQ,pct,nE,tobs,xobs,options,lhpop)
% Modified version of plotUQ_LV_timeseries for Toy Problem Data

rng(1,'twister')

iEE = 1;
% CHANGE 1: Increase threshold from 10 to 1000. 
% Your toy prey starts at 10; the original code would have discarded all your data.
che = 1000; 

if nE == 1
    for iE = 1:nUQ
        XiUQ = Xibs(:,:,iE);
        [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);
        if max(max(xSINDYb)) < che
            estimateUQ1(1:size(xSINDYb,1),iEE) = xSINDYb(:,1);
            estimateUQ2(1:size(xSINDYb,1),iEE) = xSINDYb(:,2);
            iEE = iEE + 1;
        end
    end
else
    for iE = 1:nUQ
        iMod = randi(size(Xibs,3),1,nE);
        XiUQ = mean(Xibs(:,:,iMod),3);
        [~,xSINDYb]=ode45(@(t,x)sparseGalerkin(t,x,XiUQ,polys),tspan,x0,options);
        if max(max(xSINDYb)) < che
            estimateUQ1(1:size(xSINDYb,1),iEE) = xSINDYb(:,1);
            estimateUQ2(1:size(xSINDYb,1),iEE) = xSINDYb(:,2);
            iEE = iEE + 1;
        end
    end
end

[tspanSINDyXi,xSINDYXi]=ode45(@(t,x)sparseGalerkin(t,x,Xi,polys),tspan,x0,options);

pctL = (100-pct)/2;
pctH = (100+pct)/2;
prcTT1 = prctile(estimateUQ1,[pctL 50 pctH],2);
prcTT2 = prctile(estimateUQ2,[pctL 50 pctH],2);

%% Plotting logic (Colors preserved from original)
C10 = [0 128 255]/255;
mymap3 = [255,237,160; 254,178,76; 240,59,32]./255;
cLynx = C10; cHare = mymap3(3,:);
fos = 14; fosT = 14; ms = 8; lw = 1.0; lwlb = 1.0;
cGrey = 0.6*[1 1 1]; faceAlpha = 0.2; lslb = '-';

% CHANGE 2: Remove the +1900 year offset logic
tTrue = tobs; 
tspanPlot = tspan;
tspanSINDyXiPlot = tspanSINDyXi;

figure('Position', [10 10 600 380])

% Ghost plot for legend
plot([0 1],[-1000 -1000],'-','Color',cHare,'LineWidth',lwlb); hold on
plot([0 1],[-1000 -1000],'-','Color',cLynx,'LineWidth',lwlb); hold on
plot([0 1],[-1000 -1000],'x','Color',cGrey,'LineWidth',lw,'MarkerSize',ms); hold on
plot([0 1],[-1000 -1000],lslb,'Color',cGrey,'LineWidth',lwlb); hold on
patch([0 1 1 0], [-1000 -1000 -1001 -1001], cGrey, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on

% CHANGE 3: Dynamic scaling factors based on your toy data matrix (lhpop)
n1 = std(lhpop(3,:)'); 
n2 = std(lhpop(2,:)');

% Actual data plotting
patch([tspanPlot'; flipud(tspanPlot')], [prcTT1(:,1)*n1; flipud(prcTT1(:,3)*n1)], cHare, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xobs(:,1)*n1,'x','Color',cHare,'LineWidth',lw,'MarkerSize',ms); hold on
plot(tspanSINDyXiPlot,xSINDYXi(:,1)*n1,lslb,'Color',cHare,'LineWidth',lwlb); hold on
patch([tspanPlot'; flipud(tspanPlot')], [prcTT2(:,1)*n2; flipud(prcTT2(:,3)*n2)], cLynx, 'FaceAlpha',faceAlpha, 'EdgeColor','none'); hold on
plot(tTrue,xobs(:,2)*n2,'x','Color',cLynx,'LineWidth',lw,'MarkerSize',ms); hold on
plot(tspanSINDyXiPlot,xSINDYXi(:,2)*n2,lslb,'Color',cLynx,'LineWidth',lwlb); hold on

% CHANGE 4: Update Axis limits for 0-20 time range and toy population scale
set(gca,'ticklabelinterpreter','latex','FontSize',fos)
ylabel('Population','interpreter','latex','FontSize',fos)
xlabel('Time','interpreter','latex','FontSize',fos)
% Dynamic limits:
ylim([-1 max(max(xobs*n1))*1.5]) 
xlim([min(tspanPlot)-0.5 max(tspanPlot)+0.5])

legend({'Prey (Hare)','Predator (Lynx)','Obs. data','LB-SINDy',sprintf('%d\\%% conf.',pct)},...
    'Location','NorthEast','interpreter','latex','FontSize',fos)
title('Toy Problem: Observed vs. SINDy Ensemble','interpreter','latex','Fontsize',fosT)