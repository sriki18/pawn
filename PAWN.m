function [KS,xvals,y_u, y_c, par_u, par_c, ft] = PAWN(model, p, lb, ub, ...
Nu, n, Nc, npts, seed)
%PAWN Run PAWN for Global Sensitivity Analysis of a supplied model
%   [KS,xvals,y_u, y_c, par_u, par_c, ft] = PAWN(model, p, lb, ub, ...
% Nu, n, Nc, npts, seed)
%
%     model : A function that takes a vector of parameters x and a 
%             structure p as input to provide the model output as a scalar.
%             The vector holds parameters we need sensitivity of and p
%             holds all other parameters required to run the model.
%     lb : A vector (1xM) of the lower bounds for each parameter
%     ub : A vector (1xM) of the upper bounds for each parameter
%     Nu : Number of samples of the unconditioned parameter space
%     n : Number of conditioning values in each dimension
%     Nc : Number of samples of the conditioned parameter space
%     npts : Number of points to use in kernel density estimation
%     seed : Random number seed
%
% ## Refernces
% [1]: Pianosi, F., Wagener, T., 2015. A simple and efficient method for 
% global sensitivity analysis based on cumulative distribution functions. 
% Environ. Model. Softw. 67, 1-11. doi:10.1016/j.envsoft.2015.01.004


    % Initializations
    M = length(lb); % Number of parameters
    y_u = nan(Nu,1); % Ouput of unconditioned simulations
    y_c = n * Nc * M; % Output of conditioned simulations
    KS = nan(M,n); % Kolmogorov-Smirnov statistic
    xvals = nan(M,n); % Container for conditioned samples
    ft = nan(M*n, npts); % CDF container
    
    rng(seed); % Set random seed
    % Containers for parameters
    par_u = bsxfun(@plus, lb, bsxfun(@times, rand(Nu, M), (ub-lb)));                 
    par_c = bsxfun(@plus, lb, bsxfun(@times, rand(M*Nc*n, length(lb)), (ub-lb)));
    % Create the conditioned samples
    for ind=1:M
        for ind2=1:n
            xvals(ind,ind2) = lb(ind) + rand*(ub(ind)-lb(ind));
            [(ind-1)*Nc*n+(ind2-1)*Nc+1:(ind-1)*Nc*n+ind2*Nc]
            par_c((ind-1)*Nc*n+(ind2-1)*Nc+1:...
                (ind-1)*Nc*n+ind2*Nc,ind) = xvals(ind,ind2);
        end
    end
    
    % Evaluate model output of unconditioned samples, can parallelize by
    % commenting out the parfor line and adding a comment to the for line.
%     parfor ind=1:Nu
    for ind=1:Nu
        y_u(ind) = model(par_u(ind,:), p);
    end
    
    % Evaluate model output of conditioned samples, can parallelize by
    % commenting out the parfor line and adding a comment to the for line.
%     parfor ind=1:length(par_c)
    for ind=1:length(par_c)
        y_c(ind) = model(par_c(ind,:), p);
    end
    
    % Find bounds of the model outputs
    m1 = min([y_c, y_u']);
    m2 = max([y_c, y_u']);
    % Evaluate the CDF with kernel density for unconditioned samples
    [f,~] = ksdensity(y_u, linspace(m1,m2,npts), 'Function', 'cdf');
    
    % Evaluate the CDF with kernel density for conditioned samples and use
    % that to find the KS statistic (Eqn 4 in the paper). 
    for ind=1:M
        for ind2=1:n
            % Temporarily store the current conditioned samples
            yt = y_c((ind-1)*Nc*n+(ind2-1)*Nc+1:(ind-1)*Nc*n+ind2*Nc);
            [ft((ind-1)*n+ind2,:),~] = ksdensity(yt, linspace(m1, m2,...
                                                 npts), 'Function', 'cdf');
            KS(ind,ind2) = max(abs(ft((ind-1)*n+ind2,:)-f)); % Eqn 4
        end
    end
end