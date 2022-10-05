function [ dst ] = bilateralFilter( src, k, sigma_s, sigma_r )
%bilateralFilter Applies Bilateral Filter on grayscale double image.
%   src: Source image
%   dst: Filtered image
%   k: Mask size
%   sigma_s: Closeness function parameter
%   sigma_r: Similarity function parameter

if mod(k,2)==0
    error('k should be odd');
end
if sigma_s<=0
    error('sigma_s should be positive');
end
if sigma_r<=0
    error('sigma_r should be positive');
end

dst = zeros(size(src));
% Closeness function
c = fspecial('gaussian', [k k], sigma_s);

w = floor(k/2);   % Precompute for efficiency
for i=1+w:size(src,1)-w
    for j=1+w:size(src,2)-w
        % Calculate region
        Q = src(i-w:i+w, j-w:j+w);
        
        % Similarity function
        s = exp(-(src(i,j)-Q).^2) ./ (2*sigma_r^2);
        
        kernel = c.*s;
        dst(i,j) = sum(kernel(:).*Q(:)) ./ sum(kernel(:));
    end
end

end

