clc;
close all;
clear;

videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
frameRate = vidObj.FrameRate;
numFrames = vidObj.NumFrames;
rows = vidObj.Height;
cols = vidObj.Width;
originalFrames = zeros(rows, cols, numFrames, 'uint8');
mean_Frames = zeros(rows, cols, numFrames, 'uint8');
median_Frames = zeros(rows, cols, numFrames, 'uint8');
weighted_Frames = zeros(rows, cols, numFrames, 'uint8');
diffFrames = zeros(rows, cols, numFrames, 'uint8');
for_mean_Frames = zeros(rows, cols, numFrames, 'uint8');
for_median_Frames = zeros(rows, cols, numFrames, 'uint8');
for_weighted_Frames = zeros(rows, cols, numFrames, 'uint8');
for_diff_Frames = zeros(rows, cols, numFrames, 'uint8');
for i = 1:numFrames
    originalFrames(:, :, i) = rgb2gray(read(vidObj, i));
end
Threshold = [10, 50, 100];
figure;
plotIndex = 1;
for t = 1:length(Threshold)
    T = Threshold(t);  
    mean_Frames(:, :, 1:3) = originalFrames(:, :, 1:3);
    median_Frames(:, :, 1:3) = originalFrames(:, :, 1:3);
    weighted_Frames(:, :, 1:3) = originalFrames(:, :, 1:3);
    diffFrames(:, :, 1) = originalFrames(:, :, 1);
    diffFrames(:, :, 2) = originalFrames(:, :, 2)-originalFrames(:, :, 1);
    diffFrames(:, :, 3) = originalFrames(:, :, 3)-originalFrames(:, :, 2);
    for i = 4:numFrames
      mean_Frames(:, :, i) = uint8(0.33 * originalFrames(:, :, i-1) + 0.33 * originalFrames(:, :, i-2)+ 0.33* originalFrames(:,:,i-3));
      median_Frames(:, :, i) = median(cat(3, originalFrames(:, :, max(1, i-2)), originalFrames(:, :, i-2), originalFrames(:, :, i-3)), 3);
      weighted_Frames(:, :, i) = uint8(0.5* originalFrames(:, :, i-1) + 0.3*originalFrames(:, :, i-2)+0.2*originalFrames(:, :, i-3));
      diffFrames(:, :, i) = uint8(abs(originalFrames(:, :, i) - originalFrames(:, :, i-1)));
    end
 for f = 2:numFrames
     for_mean_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(mean_Frames(:, :, f))) > T) * 255);
     for_median_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(median_Frames(:, :, f))) > T) * 255);
     for_weighted_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(weighted_Frames(:, :, f))) > T) * 255);
     for_diff_Frames(:, :, f) = uint8((diffFrames(:, :, f) > T) * 255);
 end
  subplot(4, 3, t), imshow(for_diff_Frames(:,:,165)), title(sprintf("Diff-app T=%d", T));
  subplot(4, 3, 3+t), imshow(for_mean_Frames(:,:,165)), title(sprintf("Mean-3 T=%d", T));
  subplot(4, 3, 6+t), imshow(for_median_Frames(:,:,165)), title(sprintf("Median-3 T=%d ", T));
  subplot(4, 3, 9+t), imshow(for_weighted_Frames(:,:,165)), title(sprintf("Weighted-3 T=%d ", T));
end

Threshold = {10, 50, 100};
figure();
for t = 1:3
    T = Threshold{t};
    mean_Frames(:, :, 1:5) = originalFrames(:, :, 1:5);
    median_Frames(:, :, 1:5) = originalFrames(:, :, 1:5);
    weighted_Frames(:, :, 1:5) = originalFrames(:, :, 1:5);
    for i = 6:numFrames
        mean_Frames(:, :, i) = uint8(0.2 * originalFrames(:, :, i-1) +0.2 * originalFrames(:, :, i-2) + 0.2 * originalFrames(:, :, i-3) + 0.2 * originalFrames(:, :, i-4) + 0.2 * originalFrames(:, :, i-5));
        median_Frames(:, :, i) = median(cat(3, originalFrames(:, :, i-1),originalFrames(:, :, i-2),originalFrames(:, :, i-3), originalFrames(:, :, i-4), originalFrames(:, :, i-5)), 3);
        weighted_Frames(:, :, i) = uint8(0.4 * originalFrames(:, :, i-1) + 0.3 * originalFrames(:, :, i-2) + 0.1 * originalFrames(:, :, i-3) + 0.1 * originalFrames(:, :, i-4) + 0.1 * originalFrames(:, :, i-5));
    end
    for f = 2:numFrames
        for_mean_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(mean_Frames(:, :, f))) > T) * 255);
        for_median_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(median_Frames(:, :, f))) > T) * 255);
        for_weighted_Frames(:, :, f) = uint8((abs(double(originalFrames(:, :, f)) - double(weighted_Frames(:, :, f))) > T) * 255);
    end
    subplot(3, 3, 3*t-2);  imshow(for_mean_Frames(:,:,165));    title(sprintf("Mean-5, T=%d", T));
    subplot(3, 3, 3*t-1);  imshow(for_median_Frames(:,:,165));  title(sprintf("Median-5, T=%d", T));
    subplot(3, 3, 3*t);    imshow(for_weighted_Frames(:,:,165)); title(sprintf("Weighted-5, T=%d", T));       
end
backgroundFrames2 = zeros(rows, cols, numFrames, 'uint8');
backgroundFrames5 = zeros(rows, cols, numFrames, 'uint8');
backgroundFrames8 = zeros(rows, cols, numFrames, 'uint8');
for_backgroundFrames2 = zeros(rows, cols, numFrames, 'uint8');
for_backgroundFrames5 = zeros(rows, cols, numFrames, 'uint8');
for_backgroundFrames8 = zeros(rows, cols, numFrames, 'uint8');

backgroundFrames2(:, :, 1) = originalFrames(:, :, 1);
backgroundFrames5(:, :, 1) = originalFrames(:, :, 1);
backgroundFrames8(:, :, 1) = originalFrames(:, :, 1);
for f = 2:numFrames
    backgroundFrames2(:, :, f) = uint8(0.2 * originalFrames(:, :, f-1) + 0.8* backgroundFrames2(:, :, f-1));
    backgroundFrames5(:, :, f) = uint8(0.5 * originalFrames(:, :, f-1) + 0.5* backgroundFrames5(:, :, f-1));
    backgroundFrames8(:, :, f) = uint8(0.8 * originalFrames(:, :, f-1) + 0.2* backgroundFrames8(:, :, f-1));
end
for f = 2:numFrames
    for_backgroundFrames2(:, :, f) = uint8(abs(double(originalFrames(:, :, f)) - double(backgroundFrames2(:, :, f-1))));
    for_backgroundFrames5(:, :, f) = uint8(abs(double(originalFrames(:, :, f)) - double(backgroundFrames5(:, :, f-1))));
    for_backgroundFrames8(:, :, f) = uint8(abs(double(originalFrames(:, :, f)) - double(backgroundFrames8(:, :, f-1))));
end
figure();
subplot(1, 3, 1), imshow(for_backgroundFrames2(:, :, 165)); title(sprintf("Recursive-0.2"));
subplot(1, 3, 2), imshow(for_backgroundFrames5(:, :, 165)); title(sprintf("Recursive-0.5"));
subplot(1, 3, 3), imshow(for_backgroundFrames8(:, :, 165)); title(sprintf("Recursive-0.8")); 

