function [ F_even, F_odd ] = getOrientationFilters ( theta, sigma, kernel_size, m )
% Gets the even and odd oriented energy filters.
% Based on the Gabor filter proposed in "2D feature detection via local
% energy" [Robbins and Owens].
%
% The streched Gabor function shows the radial variation, and the power of
% a cosine shows the angular variation.
%
% theta = requested orientation
% sigma = sigma of gaussian term
% kernel_size = size of the filter
% m = exponent of the cosine.
    
    if nargin < 4
        m=2;
    end

    k = floor(kernel_size/2);
    x = repmat([-k:k],[kernel_size,1]);
    y = repmat([k:-1:-k]',[1,kernel_size]);
    r = sqrt(x.^2+y.^2);
    phi = atan2(y,x);

    %% Get optimal value of lambda
    % Lambda is a parameter which makes the integral of the even-symmetric
    % orientation function equal to zero. Use the bisection method to
    % numerically calculate the value of lambda.
    
    bounds = [0 20];
    lambda = mean(bounds);
    v0 = 1;
    
    % Search for value of lambda
    while(1)
        xi = 0.5*exp(-lambda.*(r./sigma).^2) + 0.5;
        F = exp(-r.^2./(2*sigma^2)) .* cos(2*pi*v0.*r.*xi);
        integ = cumtrapz(F(:));     % Cumulative integral

        % End condition
        if abs(integ(end)) < 0.001   % Take last value in integ since it is the integral of the entire signal
            break;
        end
        
        % Bisection method to find lambda such that integral=0
        if integ(end) > 0
            bounds(1) = lambda;
            lambda = mean(bounds);
        else 
            bounds(2) = lambda;
            lambda = mean(bounds);
        end
    end
    
    % Even-Symmetric Filter
    F_even = exp(-r.^2./(2*sigma^2)) .* cos(2*pi*v0.*r.*xi) .* cos(theta - phi).^(2*m);
    % Odd-Symmetric Filter
    F_odd = exp(-r.^2./(2*sigma^2)) .* sin(2*pi*v0.*r.*xi) .* sin(theta - phi).^(2*m);

end
