clc;
close all;
clear;

videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
frameRate = vidObj.FrameRate;
numFrames = 100;
rows = vidObj.Height;
cols = vidObj.Width;

originalFrames = zeros(rows, cols, numFrames, 'uint8');
noisyFrames = zeros(rows, cols, numFrames, 'uint8');
mean_spatial_3 = zeros(rows, cols, numFrames, 'uint8');
median_spatial_3= zeros(rows, cols, numFrames, 'uint8');
weighted_spatial_3 = zeros(rows, cols, numFrames, 'uint8');
mean_temporal_3 = zeros(rows, cols, numFrames, 'uint8');
median_temporal_3= zeros(rows, cols, numFrames, 'uint8');
weighted_temporal_3 = zeros(rows, cols, numFrames, 'uint8');
mean_s_t_3 = zeros(rows, cols, numFrames, 'uint8');
median_s_t_3= zeros(rows, cols, numFrames, 'uint8');
weighted_s_t_3 = zeros(rows, cols, numFrames, 'uint8');
mean_spatial_5 = zeros(rows, cols, numFrames, 'uint8');
median_spatial_5= zeros(rows, cols, numFrames, 'uint8');
weighted_spatial_5 = zeros(rows, cols, numFrames, 'uint8');
mean_temporal_5 = zeros(rows, cols, numFrames, 'uint8');
median_temporal_5= zeros(rows, cols, numFrames, 'uint8');
weighted_temporal_5 = zeros(rows, cols, numFrames, 'uint8');
mean_s_t_5 = zeros(rows, cols, numFrames, 'uint8');
median_s_t_5= zeros(rows, cols, numFrames, 'uint8');
weighted_s_t_5 = zeros(rows, cols, numFrames, 'uint8');
tic;
for i = 1:numFrames
    originalFrames(:, :,i ) = rgb2gray(read(vidObj, i));
    noisyFrames(:, :, i) = imnoise(originalFrames(:, :, i), 'salt & pepper', 0.02);  
    %noisyFrames(:, :, i) = imnoise(originalFrames(:, :, i), 'gaussian', 0.02);
end
toc;
for i = 1: numFrames
    frame = noisyFrames(:,:,i);
    [rows, cols] = size(frame);
    weighted_kernel = [1,1,1; 1,2,1; 1,1,1] / 10;  
    for r = 2:rows-1
        for c= 2:cols-1
            window = frame(r-1:r+1, c-1:c+1);
            mean_spatial_3(r,c,i) = mean(window(:));
            median_spatial_3(r,c,i) = median(window(:));
            weighted_spatial_3(r, c, i) = sum(sum(double(window) .* double(weighted_kernel)));
        end
    end
end
for i = 1: numFrames
    frame = noisyFrames(:,:,i);
    weighted_kernel_5x5 = [1,1,1,1,1; 1,1,2,1,1; 1,2,4,2,1; 1,1,2,1,1; 1,1,1,1,1] / 24;  
    for r = 3:rows-2
        for c= 3:cols-2
            window = frame(r-2:r+2, c-2:c+2);
            mean_spatial_5(r,c,i) = mean(window(:));
            median_spatial_5(r,c,i) = median(window(:));
            weighted_spatial_5 (r, c, i) = sum(sum(double(window) .* double(weighted_kernel_5x5)));
        end
    end
end
mean_temporal_3(:,:,1)= noisyFrames(:,:,1);
median_temporal_3(:,:,1)= noisyFrames(:,:,1);
weighted_temporal_3(:,:,1)= noisyFrames(:,:,1);
for k = 2:numFrames-1
     mean_temporal_3(:,:, k) = 0.33*noisyFrames(:,:,k-1) + 0.33*noisyFrames(:,:,k)+ 0.33*noisyFrames(:,:,k+1);
     median_temporal_3(:,:, k) = median(noisyFrames(:,:, k-1:k+1), 3);
     weighted_temporal_3(:,:, k) = 0.25*noisyFrames(:,:,k-1) + 0.5*noisyFrames(:,:,k)+ 0.25*noisyFrames(:,:,k+1);

end
mean_temporal_3(:,:,numFrames)= noisyFrames(:,:,numFrames);
median_temporal_3(:,:,numFrames)= noisyFrames(:,:,numFrames);
weighted_temporal_3(:,:,numFrames)= noisyFrames(:,:,numFrames);
mean_temporal_5(:,:,1)= noisyFrames(:,:,1);median_temporal_5(:,:,1)= noisyFrames(:,:,1);weighted_temporal_5(:,:,1)= noisyFrames(:,:,1);
mean_temporal_5(:,:,1)= noisyFrames(:,:,1);median_temporal_5(:,:,2)= noisyFrames(:,:,2);weighted_temporal_5(:,:,2)= noisyFrames(:,:,1);
for k = 3:numFrames-2
     mean_temporal_5(:,:, k) = 0.2*noisyFrames(:,:,k-2) + 0.2*noisyFrames(:,:,k-1)+ 0.2*noisyFrames(:,:,k)+0.2*noisyFrames(:,:,k+1)+ 0.2*noisyFrames(:,:,k+2);
     median_temporal_5(:,:, k) = median(noisyFrames(:,:, k-2:k+2), 3);
     weighted_temporal_5(:,:, k) = 0.1*noisyFrames(:,:,k-2) + 0.2*noisyFrames(:,:,k-1)+ 0.4*noisyFrames(:,:,k)+0.2*noisyFrames(:,:,k+1)+0.1*noisyFrames(:,:,k+2);
end
mean_temporal_5(:,:,numFrames-1)= noisyFrames(:,:,numFrames-1);median_temporal_5(:,:,numFrames-1)= noisyFrames(:,:,numFrames-1);weighted_temporal_5(:,:,numFrames-1)= noisyFrames(:,:,numFrames-1);
mean_temporal_5(:,:,numFrames)= noisyFrames(:,:,numFrames);median_temporal_5(:,:,numFrames)= noisyFrames(:,:,numFrames);weighted_temporal_5(:,:,numFrames)= noisyFrames(:,:,numFrames);

for i =  2: numFrames-1
   for r = 2 : rows-1 
       for c = 2: cols-1 
          median_s_t_3(r, c, i) = median(median(median(cat(3, ...
          noisyFrames(r-1:r+1, c-1:c+1, i-1), ...
          noisyFrames(r-1:r+1, c-1:c+1, i), ...
          noisyFrames(r-1:r+1, c-1:c+1, i+1)), 3)));
       end 
   end
end
for i = 2: numFrames-1
     weighted_s_t_3(:,:,i) = 0.25* noisyFrames(:,:,i-1)+ 0.5* weighted_spatial_3(:,:,i)+0.25*noisyFrames(:,:,i+1);
      mean_s_t_3(:,:,i) = 0.33*mean_spatial_3(:,:,i-1) + 0.33*mean_spatial_3(:,:,i)+ 0.33*mean_spatial_3(:,:,i+1);
end 
for i =  3: numFrames-2
    mean_s_t_5(:,:,i) = 0.2*mean_spatial_5(:,:,i-2)+0.2*mean_spatial_5(:,:,i-1) + 0.2*mean_spatial_5(:,:,i)+ 0.2*mean_spatial_5(:,:,i+1)+0.2*mean_spatial_5(:,:,i-2);
    weighted_s_t_5(:,:,i) = 0.1* noisyFrames(:,:,i-2)+ 0.2* weighted_spatial_3(:,:,i-1)+0.4* weighted_spatial_5(:,:,i)+ 0.2* weighted_spatial_3(:,:,i+1)+0.1*noisyFrames(:,:,i+2);
end
for i =  3: numFrames-2
   for r = 3 : rows-2
       for c = 3: cols-2 
        median_s_t_5(r, c, i) = median(median(median(cat(3, ...
        noisyFrames(r-2:r+2, c-2:c+2, i-2), ...
        noisyFrames(r-2:r+2, c-2:c+2, i-1), ...
        noisyFrames(r-2:r+2, c-2:c+2, i), ...
        noisyFrames(r-2:r+2, c-2:c+2, i+1), ...
        noisyFrames(r-2:r+2, c-2:c+2, i+2)), 3), 'all'));
       end 
   end
end
filter_list_3 = {mean_spatial_3,median_spatial_3, weighted_spatial_3,mean_temporal_3,median_temporal_3, weighted_temporal_3,mean_s_t_3,median_s_t_3, weighted_s_t_3};
filter_list_5 = {mean_spatial_5, median_spatial_5, weighted_spatial_5,mean_temporal_5, median_temporal_5, weighted_temporal_5,mean_s_t_5, median_s_t_5, weighted_s_t_5};
mse_values_3 = zeros(numFrames, 9);
psnr_values_3 = zeros(numFrames, 9);
mse_values_5 = zeros(numFrames, 9);
psnr_values_5 = zeros(numFrames, 9);

for i = 1:9
    for k = 1:numFrames
        mse_values_3(k,i) = mean((double(originalFrames(:,:,k)) - double(filter_list_3{i}(:,:,k))).^2, 'all');
        mse_values_5(k,i) = mean((double(originalFrames(:,:,k)) - double(filter_list_5{i}(:,:,k))).^2, 'all');
        if mse_values_3(k,i) == 0
            psnr_values_3(k,i) = Inf; % Avoid log(0)
        else
            psnr_values_3(k,i) = 10 * log10(255^2 / mse_values_3(k,i));
        end
        if mse_values_5(k,i) == 0
            psnr_values_5(k,i) = Inf;
        else
            psnr_values_5(k,i) = 10 * log10(255^2 / mse_values_5(k,i));
        end
    end
end

filter_name_3 = {'Mean Spatial 3', 'Median Spatial 3', 'Weighted Spatial 3', ...
                  'Mean Temporal 3', 'Median Temporal 3', 'Weighted Temporal 3', ...
                  'Mean Spatial-Temporal 3', 'Median Spatial-Temporal 3', 'Weighted Spatial-Temporal 3'};
filter_name_5 = {'mean_spatial_5','median_spatial_5', 'weighted_spatial_5','mean_temporal_5', 'median_temporal_5', 'weighted_temporal_5','mean_Spatial-Temporal_5','median_Spatial-Temporal_5','weighted_Spatial-Temporal_5'};

figure();
subplot(1,2,1); imshow(originalFrames(:,:,5)); title("original");
subplot(1,2,2); imshow(noisyFrames(:,:,5)); title("salt-pepper-nosied");
figure();
for i = 1: 5
    avg_mse_3 = mean(mse_values_3(:,i));
    avg_mse_5 = mean(mse_values_5(:,i));
    avg_psnr_3 = 10 * log10(255^2 / avg_mse_3);
    avg_psnr_5 =  10 * log10(255^2 / avg_mse_5);
    subplot(5,2,2*i-1); plot(1:numFrames, psnr_values_3(:,i), '-');
    title(sprintf('%s\nMSE: %f, PSNR: %f dB', strrep(filter_name_3{i}, '_', ' '),avg_mse_3 , avg_psnr_3), 'FontSize', 10);
    subplot(5,2,2*i); plot(1:numFrames, psnr_values_5(:,i), '-');
    title(sprintf('%s\nMSE: %f, PSNR: %f dB', strrep(filter_name_5{i}, '_', ' '),avg_mse_5 , avg_psnr_5), 'FontSize', 10);
end
figure()
for i = 6 :9
    avg_mse_3 = mean(mse_values_3(:,i));
    avg_mse_5 = mean(mse_values_5(:,i));
    avg_psnr_3 = 10 * log10(255^2 / avg_mse_3);
    avg_psnr_5 =  10 * log10(255^2 / avg_mse_5);
    subplot(4,2,2*i-11); plot(1:numFrames, psnr_values_3(:,i), '-');
    title(sprintf('%s\nMSE: %f, PSNR: %f dB', strrep(filter_name_3{i}, '_', ' '),avg_mse_3 , avg_psnr_3), 'FontSize', 10);
    subplot(4,2,2*i-10); plot(1:numFrames, psnr_values_5(:,i), '-');
    title(sprintf('%s\nMSE: %f, PSNR: %f dB', strrep(filter_name_5{i}, '_', ' '),avg_mse_5 , avg_psnr_5), 'FontSize', 10);
end
figure();
for f = 1:numel(filter_list_3)
    subplot(3, 3, f);
    imshow(filter_list_3{f}(:,:,5));   % Convert string to variable
    title(strrep(filter_name_3{f}, '_', ' '), 'FontSize', 10);  % Replace underscores with spaces for better readability
end
figure();
for f = 1:numel(filter_list_5)
    subplot(3, 3, f);
    imshow(filter_list_5{f}(:,:,5));   % Convert string to variable
    title(strrep(filter_name_5{f}, '_', ' '), 'FontSize', 10); % Replace underscores with spaces for better readability
end
