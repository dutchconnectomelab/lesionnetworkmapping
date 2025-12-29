function [h, crit_p, adj_p] = fdr(pvals, q)
    % FDR Benjamini–Hochberg false discovery rate (BH-FDR) correction.
    %
    %   [h, crit_p, adj_p] = fdr(pvals)
    %   [h, crit_p, adj_p] = fdr(pvals, q)
    %
    %   Applies the Benjamini–Hochberg procedure to control the expected false
    %   discovery rate at level q. Works on vectors or matrices; computations are
    %   performed on the flattened p-value list and then reshaped back.
    %
    % Inputs:
    %   pvals   : Vector or matrix of raw p-values.
    %   q       : FDR rate (scalar), default 0.05.
    %
    % Outputs:
    %   h       : Logical array, same size as pvals (true = significant under BH at q).
    %   crit_p  : Critical p-value threshold (largest p passing BH); 0 if none pass.
    %   adj_p   : BH-adjusted p-values (q-values), same size as pvals.
    %
    % Example
    %   [h, crit_p, qvals] = fdr(pvals, 0.05);

    if nargin < 2
        q = 0.05;
    end

    % Flatten to vector
    p = pvals(:);
    m = length(p);  % total number of tests

    % Sort p-values and get their original indices
    [sorted_p, sort_idx] = sort(p);
    unsort_idx = zeros(size(sort_idx));
    unsort_idx(sort_idx) = 1:m;

    % Compute FDR thresholds
    bh_thresholds = (1:m)' / m * q;

    % Find largest p-value that is below the threshold
    below_threshold = sorted_p <= bh_thresholds;
    max_id = find(below_threshold, 1, 'last');

    if isempty(max_id)
        crit_p = 0;
        h = false(size(p));
    else
        crit_p = sorted_p(max_id);
        h = p <= crit_p;
    end

    % Compute adjusted p-values
    adj_p = zeros(size(p));
    min_adj_p = 1;
    for i = m:-1:1
        min_adj_p = min(min_adj_p, sorted_p(i) * m / i);
        adj_p(sort_idx(i)) = min_adj_p;
    end

    % Reshape to original shape
    h = reshape(h, size(pvals));
    adj_p = reshape(adj_p, size(pvals));
end