%% Trilateral Filtering for Biomedical Images
% By Claudio Alkazzi: 2096201
%
% Summary:
% This program takes an image of a brain containing a tumor and corrupts it
% with different levels of noise. It then performs bilateral and trilateral
% filtering on all the noisy images and plots the MSE for each filter and
% for each level of noise.
% 
% Outputs:
%   - A plot of the MSE
%   - The resulting images in a folder called data. This folder already
%   contains the output of the code. Running this code again will override
%   the existing files.
%
% User-defined Functions:
%   - bilateralFilter.m: Perfors bilateral filtering on an input image.
%   - trilateralFilter.m: Performs trilateral filtering on an input image.
%   - getOrientationFilters.m: Helper function for
%   getLocalOrientationEnergy.m
%   - getLocalOrientationEnergy: Helper function for trilateralFilter.m
%
%% Preprocess Input Image
src = imread('../data/brain_tumor.png');
if size(src,3) == 3
    src = rgb2gray(src);
end
src = im2double(src);
img = src;  % Image without noise
% Variance of the Gaussian noise
noise = [0.01, 0.02, 0.03, 0.05, 0.1];
mse = zeros(3,length(noise));

%% Perform Bilateral and Trilateral Filtering
for count = 1:length(noise)
    % Corrupt image with noise and output to jpg image
    src = imnoise(img,'gaussian',0,noise(count));
    filename = strcat('src_',num2str(noise(count)),'.jpg');
    filepath = strcat('../data/output/',filename);
    imwrite(src,filepath);
    
    % Bilateral Filter
    dst = src;
    for i=1:3
        dst = bilateralFilter(dst, 5, 16, 0.2);
    end
    filename = strcat('bf_',num2str(noise(count)),'.jpg');
    filepath = strcat('../data/output/',filename);
    imwrite(dst,filepath);

    % Trilateral Filter
    dst1 = src;
    for i=1:3
        dst1 = trilateralFilter(dst1, 3, 8, 0.2, 0.2);
    end
    filename = strcat('tf_',num2str(noise(count)),'.jpg');
    filepath = strcat('../data/output/',filename);
    imwrite(dst1,filepath);
    
    % MSE Calculation
    % Store MSE for noisy image, bilateral filter, and trilateral filter
    % compared to original image.
    mse(:,count) = [immse(src,img) immse(dst,img) immse(dst1,img)];
end

%% Plot MSE Curves
figure, hold on
for i = 1:3
    plot(noise,mse(i,:));
end
title('MSE');
xlabel('Variance of Gaussian Noise');
ylabel('MSE value');
legend('Noisy Image', 'BF', 'TF');
hold off


