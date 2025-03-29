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
figure();
subplot(1,2,1), imshow(gray1); title('Frame 1');
subplot(1,2,2), imshow(gray2); title('Frame 2');
tic;
for i = 1:blocksize:m-blocksize
    for j = 1:blocksize:n-blocksize
        block = gray1(i:i+blocksize-1, j:j+blocksize-1);
        minMAD = Inf;  best_x = i; best_y = j;
        upper_row = max(1, i - search_range);
        lower_row = min(i + search_range, m - blocksize + 1);
        left_col = max(1, j - search_range);
        right_col = min(j + search_range, n - blocksize + 1);
        for x = upper_row:lower_row
            for y = left_col:right_col
                search_block = gray2(x:x+blocksize-1, y:y+blocksize-1);
                MAD = sum(abs(double(block) - double(search_block)), 'all');
                if MAD < minMAD
                    minMAD = MAD; best_x = x;    best_y = y;
                end
            end
        end
        row_idx = floor(i/blocksize) + 1; col_idx = floor(j/blocksize) + 1;
        motionVectorsX(row_idx, col_idx) = best_y - j;
        motionVectorsY(row_idx, col_idx) = best_x - i;
        startX(row_idx, col_idx) = best_y;
        startY(row_idx, col_idx) = best_x;
    end
end
toc;
figure;
imshow(gray1);
title('Motion Vector Field on Frame 2');
hold on;
quiver(startX, startY, motionVectorsX, motionVectorsY, 'r', 'LineWidth', 1.5);
hold off;
