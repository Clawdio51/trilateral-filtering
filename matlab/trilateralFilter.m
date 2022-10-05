function [ dst ] = trilateralFilter( src, kernel_size, sigma_s, sigma_r, sigma_d )
%trilateralFilter applies Trilateral Filter on grayscale double image
%   src: Source image
%   dst: Filtered image
%   k: Mask size
%   sigma_s: Closeness function parameter
%   sigma_r: Similarity function parameter
%   sigma_d: Orientation function parameter

    dst = zeros(size(src));
    
    % sigma is a tweakable parameter for orientation filters
    sigma = sqrt(2);
    [E,phi] = getLocalOrientationEnergy(src, sigma, kernel_size);

    % Closeness function
    c = fspecial('gaussian', [kernel_size kernel_size], sigma_s);

    k = floor(kernel_size/2);   % Precompute for efficiency
    for i=1+k:size(src,1)-k
        for j=1+k:size(src,2)-k
            % Calculate region
            Q = src(i-k:i+k, j-k:j+k);

            % Similarity function
            s = exp(-(src(i,j)-Q).^2) ./ (2*sigma_r^2);
            % Orientation function
            E_p = E(i,j);
            E_q = E(i-k:i+k, j-k:j+k);
            % This approximation works since we are using 3x3 kernel. If a
            % larger kernel is needed, refer to equation 15 in "Fast Trilateral
            % Filtering" [Robbins and Owens]
            deltaE_pq = (E_p + E_q) .* (phi(i,j) ~= phi(i-k:i+k, j-k:j+k));
            d = exp(-(deltaE_pq).^2 ./ (2*sigma_d^2));

            kernel = c.*s + sqrt(c.*s.*d);
            %kernel = c.*s.*d;
            dst(i,j) = sum(kernel(:).*Q(:)) ./ sum(kernel(:));
        end
    end

end

