function r2 = compute_gp_r2(gprMdl, X, y)
% Compute R² for GP fit quality assessment
    y_pred = predict(gprMdl, X);
    ss_res = sum((y - y_pred).^2);
    ss_tot = sum((y - mean(y)).^2);
    r2 = 1 - ss_res / ss_tot;
end

