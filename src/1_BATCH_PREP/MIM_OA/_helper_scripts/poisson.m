function [y] = poisson(data_in,lambda)
    y = zeros(size(data_in));
    for i = 1:length(data_in)
        y(i) = ((lambda^data_in(i))/factorial(data_in(i)))*exp(-lambda);
    end
end
