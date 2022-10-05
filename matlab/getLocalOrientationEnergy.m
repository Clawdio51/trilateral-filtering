function [ E, phi ] = getLocalOrientationEnergy( src, sigma_d, kernel_size, m )
% Gets the local orientation energy given orientation functions
%   
% E = Local Energy
% phi = phase of a pixel as defined in Fast Trilateral Filtering [Robbins
% and Owens] (Page 5 equation 12)
%
% src = source image
% F_even = even-symmetric orientation function
% F_odd = odd-symmetric orientation function
    
    if nargin < 4
        m=2;
    end
    
    E = zeros(size(src));   % Local Energy
    k = floor(kernel_size/2);
    x = repmat([-k:k],[kernel_size,1]);
    y = repmat([k:-1:-k]',[1,kernel_size]);
    theta = atan2(y,x);
    
    % Find local energy
    theta_max = 0; % Store theta which maximizes local energy.
    max = 0; % Store maximum response
    for i = 1:size(theta,1)
        for j = 1:size(theta,2)
            t = theta(i,j);
            [F_even, F_odd] = getOrientationFilters(t,sigma_d,kernel_size,m);
            response = sqrt( (conv2(src,F_even,'same')).^2 + (conv2(src,F_odd,'same')).^2 );
            E = E + response;
            response = sum(response(:));
            
            % Get the maximum response and the value of theta which
            % maximizes it. This will be used when computing phi.
            if response > max
                max = response;
                theta_max = t;
            end
        end
    end
    
    % Used to get boundary strength between two neighboring pixels
    [F_even, F_odd] = getOrientationFilters(theta_max,2,kernel_size);
    conv2(src,F_even);
    phi = sign(conv2(src,F_even,'same'));
    phi(phi==0) = -1;
end

