clc; clear; close all;
image = imread("peacock.jpg");
I2 = rgb2ycbcr(image);
Y = double(I2(:,:,1));
Cb = double(I2(:,:,2));
Cr = double(I2(:,:,3));
subCb = imresize(Cb, 0.5, 'bilinear');
subCr = imresize(Cr, 0.5, 'bilinear');
cb_u = imresize(subCb, 2, 'bilinear');
cb_r= imresize(subCr, 2, 'bilinear');
down_image = ycbcr2rgb(uint8(cat(3,Y,cb_u,cb_r)));
figure();
subplot(2,3,1); imshow(image); title("original Image");
subplot(2,3,2); imshow(down_image); title("down_sampled");
Y = Y -128 ; subCr = subCr -128 ; subCb = subCb -128;
Qy= [16 11 10 16 24 40 51 61; 
             12 12 14 19 26 58 60 55;
             14 13 16 24 40 57 69 56;
             14 17 22 29 51 87 80 62;
             18 22 37 56 68 109 103 77;
             24 35 55 64 81 104 113 92;
             49 64 78 87 103 121 120 101;
             72 92 95 98 112 100 103 99];

Qc = [17 18 24 47 99 99 99 99;
               18 21 26 66 99 99 99 99;
               24 26 56 99 99 99 99 99;
               47 66 99 99 99 99 99 99;
               99 99 99 99 99 99 99 99;
               99 99 99 99 99 99 99 99;
               99 99 99 99 99 99 99 99;
               99 99 99 99 99 99 99 99];
[m, n] = size(Y);
[mc, nc] = size(subCb);
block_size = 8;
rows = block_size * floor(m / block_size);
cols = block_size * floor(n / block_size);
rowsc = block_size * floor(mc / block_size);
colsc = block_size * floor(nc / block_size);
qf_values = [10, 30, 60 , 90 ];

for idx = 1:4 
    QF = qf_values(idx);
    sf = (QF >= 50) * ((100 - QF) / 50) + (QF < 50) * (50 / QF);
    qsy =max(1,round(sf * Qy)) ;
    qscbcr =max(1,round(sf * Qc)) ;
    Y = Y(1:rows, 1:cols);
    subCb = subCb(1:rowsc, 1:colsc);
    subCr = subCr(1:rowsc, 1:colsc);
    blocks_Y = mat2cell(Y, repmat(block_size, 1, rows / block_size), repmat(block_size, 1, cols / block_size));
    blocks_Cb = mat2cell(subCb, repmat(block_size, 1, rowsc / block_size), repmat(block_size, 1, colsc / block_size));
    blocks_Cr = mat2cell(subCr, repmat(block_size, 1, rowsc / block_size), repmat(block_size, 1, colsc / block_size));

    for i = 1:size(blocks_Y, 1)
        for j = 1:size(blocks_Y, 2)
            block_Y = blocks_Y{i, j};
            gY{i, j} = dct2(block_Y);
            bY{i, j} = round(gY{i, j} ./ qsy);
            gcapY{i, j} = bY{i, j} .* qsy;
            g_capY = gcapY{i, j};
            fcapY{i, j} = idct2(g_capY) + 128;
        end
    end

    for i = 1:size(blocks_Cb, 1)
        for j = 1:size(blocks_Cb, 2)
            block_Cb = blocks_Cb{i, j};
            gCb{i, j} = dct2(block_Cb);
            bCb{i, j} = round(gCb{i, j} ./ qscbcr);
            gcapCb{i, j} = bCb{i, j} .* qscbcr;
            g_capCb = gcapCb{i, j};
            fcapCb{i, j} = idct2(g_capCb) + 128;
        end
    end

    for i = 1:size(blocks_Cr, 1)
        for j = 1:size(blocks_Cr, 2)
            block_Cr = blocks_Cr{i, j};
            gCr{i, j} = dct2(block_Cr);
            bCr{i, j} = round(gCr{i, j} ./ qscbcr);
            gcapCr{i, j} = bCr{i, j} .* qscbcr;
            g_capCr = gcapCr{i, j};
            fcapCr{i, j} = idct2(g_capCr) + 128;
        end
    end
    
    Y1 = cell2mat(fcapY);
    Cb1 = cell2mat(fcapCb);
    Cr1 = cell2mat(fcapCr);
    
    Cbre = imresize(Cb1, [size(Y1, 1), size(Y1, 2)], 'bilinear');
    Crre = imresize(Cr1, [size(Y1, 1), size(Y1, 2)], 'bilinear');
    
    I3 = cat(3, Y1, Cbre, Crre);
    I3 = ycbcr2rgb(uint8(I3));   
    original_size = numel(image); 
    bY_mat = cell2mat(bY);
    bCb_mat = cell2mat(bCb);
    bCr_mat = cell2mat(bCr);
    quantized_size = nnz(bY_mat) + nnz(bCb_mat) + nnz(bCr_mat);
    compression_rate = original_size / quantized_size;
    mse = mean((image(:) - I3(:)).^2);
    psnr_value = psnr(I3, image);
    subplot(2,3,idx+2), imshow(I3);
    title(sprintf('QF: %d ,MSE: %.4f\nPSNR: %.4f dB\nCompression Ratio: %.4f',QF, mse, psnr_value, compression_rate));
end
