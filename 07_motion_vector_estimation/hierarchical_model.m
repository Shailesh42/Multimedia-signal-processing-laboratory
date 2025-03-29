clc; clear; close all;
videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
gray1 = rgb2gray(read(vidObj, 165));
gray2 = rgb2gray(read(vidObj, 166));

blocksize = 16;
search_range = 16; 
levels = [4, 2, 1]; 
[m, n] = size(gray1);
figure;
for idx = 1:length(levels)
    tic;
    level = levels(idx);
    img1 = imresize(gray1, 1/level);
    img2 = imresize(gray2, 1/level); 
    [m_ds, n_ds] = size(img1);
    bs = blocksize / level; 
    sr = search_range / level; 
    motionX = zeros(floor(m_ds/bs), floor(n_ds/bs));
    motionY = zeros(floor(m_ds/bs), floor(n_ds/bs));
    for i = 1:bs:m_ds-bs
        for j = 1:bs:n_ds-bs
            block = img1(i:i+bs-1, j:j+bs-1);
            minMAD = Inf;
            best_x = i;
            best_y = j;
            for x = max(1, i-sr) : min(i+sr, m_ds-bs+1)
                for y = max(1, j-sr) : min(j+sr, n_ds-bs+1)
                    search_block = img2(x:x+bs-1, y:y+bs-1);
                    MAD = sum(abs(double(block) - double(search_block)), 'all');   
                    if MAD < minMAD
                        minMAD = MAD;
                        best_x = x;
                        best_y = y;
                    end
                end
            end
            
            row_idx = floor(i/bs) + 1;
            col_idx = floor(j/bs) + 1;  
            motionX(row_idx, col_idx) = best_y - j;
            motionY(row_idx, col_idx) = best_x - i;
        end
    end
    startX = zeros(size(motionX));
    startY = zeros(size(motionY));
    for i = 1:size(motionX,1)
        for j = 1:size(motionX,2)
            startX(i,j) = (j-1) * bs + bs / 2;
            startY(i,j) = (i-1) * bs + bs / 2;
        end
    end
    subplot(1, length(levels), idx);
    imshow(img2);
    title(['Motion Vectors (Downsampling Factor: ', num2str(level), ')']);
    hold on;
    quiver(startX, startY, motionX, motionY, 'r', 'LineWidth', 1.5);
    hold off;
    toc;
end

motionVectorsX = imresize(motionX, [floor(m/blocksize), floor(n/blocksize)]) * 2;
motionVectorsY = imresize(motionY, [floor(m/blocksize), floor(n/blocksize)]) * 2;
for i = 1:size(motionVectorsX,1)
    for j = 1:size(motionVectorsX,2)
        startX(i,j) = (j-1) * blocksize + blocksize / 2;
        startY(i,j) = (i-1) * blocksize + blocksize / 2;
    end
end

figure;
imshow(gray1);
title('Hierarchical Motion Vector Field (Full Resolution)');
hold on;
quiver(startX, startY, motionVectorsX, motionVectorsY, 'r', 'LineWidth', 1.5);
hold off;
