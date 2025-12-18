function plotUQ_LV_toy(XiE,XiTrue,Xi,lib)
% Rewritten for toy problem parameter scaling

C10 = [0 128 255]/255;
CHare = [240,59,32]./255;
nP1 = size(XiE,1); % Library terms (10)
nP2 = size(XiE,2); % States (2)

figure('Position', [50 50 500 800])
nn = 1;

for ip = 1:nP1
    for jp = 1:nP2
        xPtest = squeeze(XiE(ip,jp,:));
        
        subplot(nP1, nP2, nn)
        
        % Select color based on Prey or Predator column
        if jp == 1; Clb = CHare; else; Clb = C10; end
        
        if sum(xPtest ~= 0) > 0
            % Plot the distribution of discovered coefficients
            [f, xi_pts] = ksdensity(xPtest(xPtest~=0));
            fill(xi_pts, f, Clb, 'FaceAlpha', 0.4, 'EdgeColor', Clb); hold on
            
            % Draw the "Mean" discovered coefficient line
            mY = max(f);
            plot([Xi(ip,jp) Xi(ip,jp)], [0 mY], 'k-', 'LineWidth', 1.5);
            
            % Highlight if this term is "True" in your ODE
            if XiTrue(ip,jp) ~= 0
                patch([min(xi_pts) max(xi_pts)], [0 mY], Clb, 'FaceAlpha', 0.05);
            end
        else
            text(0.5, 0.5, '0', 'HorizontalAlignment', 'center');
        end
        
        % Dynamic X-axis: Focus on the discovered values instead of fixed [-2,2]
        if any(xPtest ~= 0)
            xlim([min(xPtest)-0.1, max(xPtest)+0.1]);
        end
        
        yticks([]);
        if jp == 1
            ylabel(lib{ip}, 'interpreter', 'latex', 'Rotation', 0, 'HorizontalAlignment', 'right');
        end
        if ip == nP1
            xlabel('Value');
        end
        nn = nn + 1;
    end
end
sgtitle('Coefficient Uncertainty Distributions');