% Input:
%     alpha        a scalar. The concentration parameter of the DP.
% Output:
%     z            a column vector. Its length is the number of the
%                  observations. The values are the indicators of the
%                  observations.
function z = dp_kmeans(data, alpha, actN, maxIter)
if nargin < 2
    alpha = 1;
end
if nargin < 3
    actN = 100;
end
if nargin < 4
    maxIter = 100;
end
%--------------------------------------------------------------------------
% STEP 1: Init
%--------------------------------------------------------------------------
z = ones(size(data, 1), 1);
centers = zeros(actN, 2);
covariances = zeros(2, 2, actN);
data_mu = mean(data);
data_std = mean(std(data));

%--------------------------------------------------------------------------
% STEP 2: Gibbs sampling
%--------------------------------------------------------------------------
for iter = 1:maxIter
    % sampling G0
    counts = histcounts(z, 1:actN+1);
    a = counts + 1;
    b = [cumsum(counts(2:end), 'reverse'), 0];
    b = b + alpha;
    V = betarnd(a, b);
    G0 = V;
    V = cumprod(1-V);
    G0(2:end) = G0(2:end) .* V(1:end-1);
    
    % sample positions

    ix = accumarray(z, 1:length(z), [], @(x){x});
    for i = 1:length(ix)
        if ~isempty(ix{i})
            dd = data(ix{i}, :);
            centers(i,:) = mean(dd);
            diff = dd - repmat(centers(i,:), length(ix{i}), 1);
            if length(ix{i}) > 3
                covariances(:, :, i) = diff' * diff / length(ix{i});
            else
                [~,covariances(:, :, i)] = sample_from_prior(data_mu, data_std);
            end
            
        else
            [centers(i,:), covariances(:,:,i)] = sample_from_prior(data_mu,...
                                                                   data_std);
        end
    end
    if length(ix) < actN
        for i = (length(ix) +1) : actN
            [centers(i,:), covariances(:,:,i)] = sample_from_prior(data_mu,...
                                                                   data_std);
        end
    end
    
    % sample z
    for i = 1:size(data, 1)
        likelihood = mvnpdf(centers, data(i,:), covariances);
        post = likelihood' .* G0;
        post = post / sum(post);
        [~, ~, z(i)] = histcounts(rand(1), [0, cumsum(post)]);
    end
    fprintf('iteration %d done \n', iter)
        
end

end

function [center, covariance] = sample_from_prior(mu, std)
center = mvnrnd(mu, std * eye(2));
covariance = eye(2) * std / 10;
end
    