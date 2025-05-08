function [model,p_value,fit_boot] = fitScatter(X,Y,options)

arguments
    X double
    Y double

    options.type string = 'original'
    options.nboot double = 1000

    options.weights double
end

X = X(:);  
Y = Y(:);

mask = ~isnan(X) & ~isnan(Y); 
X = X(mask); 
Y = Y(mask);

if strcmpi(options.type,'weighted')
    A = [X(:) ones(length(X), 1)];
    weights = options.weights(mask);
    model = lscov(A, Y, weights);
    slope = model(1);  % p(1) is the slope, p(2) is the intercept.
    fit_boot = bootstrp(options.nboot, @(idx) lscov(A, Y(idx), weights(idx)), 1:length(Y));
    p_value = sum(abs(fit_boot(:,1)) >= abs(slope)) / options.nboot;
else
    model = polyfit(X, Y, 1); slope = model(1);
    fit_boot = bootstrp(options.nboot, @(i) polyfit(X, Y(i), 1), 1:length(Y));
    p_value = sum(abs(fit_boot(:,1)) >= abs(slope))/options.nboot;
    return
end

end
