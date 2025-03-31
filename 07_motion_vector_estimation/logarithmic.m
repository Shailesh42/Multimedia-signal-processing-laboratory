clc; clear; close all;
videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
gray1 = rgb2gray(read(vidObj, 165));
gray2 = rgb2gray(read(vidObj, 166));
blocksize = 16;
[m, n] = size(gray1);
motionVectorsX = zeros(floor(m/blocksize), floor(n/blocksize));
motionVectorsY = zeros(floor(m/blocksize), floor(n/blocksize));
startX = zeros(floor(m/blocksize), floor(n/blocksize));
startY = zeros(floor(m/blocksize), floor(n/blocksize));
search_range = 16; 
tic;
for i = 1:blocksize:m-blocksize
    for j = 1:blocksize:n-blocksize
        block = gray1(i:i+blocksize-1, j:j+blocksize-1);
        minMAD = Inf;
        best_x = i;
        best_y = j;
        stepSize = search_range / 2; 
        while stepSize >= 1
            candidates = [
                best_x, best_y; best_x - stepSize, best_y; 
                best_x + stepSize, best_y; best_x, best_y - stepSize; 
                best_x, best_y + stepSize; 
                best_x - stepSize, best_y - stepSize; 
                best_x - stepSize, best_y + stepSize; 
                best_x + stepSize, best_y - stepSize; 
                best_x + stepSize, best_y + stepSize; 
            ];
            for k = 1:size(candidates, 1)
                x = candidates(k, 1);
                y = candidates(k, 2);
                if x >= 1 && x <= m - blocksize + 1 && y >= 1 && y <= n - blocksize + 1
                    search_block = gray2(x:x+blocksize-1, y:y+blocksize-1);
                    MAD = sum(abs(double(block) - double(search_block)), 'all');

                    if MAD < minMAD
                        minMAD = MAD;
                        best_x = x;
                        best_y = y;
                    end
                end
            end   
            stepSize = stepSize / 2; 
        end
        row_idx = floor(i/blocksize) + 1;
        col_idx = floor(j/blocksize) + 1;
        motionVectorsX(row_idx, col_idx) = best_y - j;
        motionVectorsY(row_idx, col_idx) = best_x - i;
        startX(row_idx, col_idx) = best_y;
        startY(row_idx, col_idx) = best_x;
    end
end
toc;
figure;
imshow(gray1);
title('Logarithmic Search');
hold on;
quiver(startX, startY, motionVectorsX, motionVectorsY, 'r', 'LineWidth', 1.5);
hold off;
