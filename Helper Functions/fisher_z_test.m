function [Z_diff, p_value] = fisher_z_test(r1, r2, n1, n2)
    % Fisher Z-test for comparing two correlation coefficients
    % Inputs:
    %   r1 - correlation coefficient for first pair (e.g., x(t) vs. y(t))
    %   r2 - correlation coefficient for second pair (e.g., x(t) vs. y(t+1))
    %   n1 - sample size for first pair
    %   n2 - sample size for second pair
    % Outputs:
    %   Z_diff - Z-score for the difference between the correlations
    %   p_value - p-value of the test

    % Fisher Z-transformation for both correlations
    Z1 = 0.5 * log((1 + r1) / (1 - r1));
    Z2 = 0.5 * log((1 + r2) / (1 - r2));
    
    % Standard error of the difference
    SE = sqrt(1 / (n1 - 3) + 1 / (n2 - 3));

    % Z-score for the difference between correlations
    Z_diff = (Z1 - Z2) / SE;
    
    % Two-tailed p-value based on the Z-distribution
    p_value = 2 * (1 - normcdf(abs(Z_diff)));

end
