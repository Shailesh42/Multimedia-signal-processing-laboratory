videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
numFrames = floor(vidObj.Duration * vidObj.FrameRate);
rows = vidObj.Height;
cols = vidObj.Width;

% Quantization matrices
q_intra = [8 16 19 22 26 27 29 34;
           16 16 22 24 27 29 34 37;
           19 22 26 27 29 34 34 38;
           22 22 26 27 29 34 37 40;
           22 26 27 29 32 35 40 48;
           26 27 29 32 35 40 48 58;
           26 27 29 34 38 46 56 89;
           27 29 35 38 46 56 69 83];

q_inter = 16 * ones(8, 8);

% Parameters
blockSize = 16;
searchRange = 4;
skipThreshold = 5;

% Frame types
frameTypes = repmat({'I','B','B','P','B','B','P','B','B'}, 1, ceil(numFrames/9));
frameTypes = frameTypes(1:numFrames);

% Preallocation
encodedFrames = cell(1, numFrames);
decodedRGBFrames = zeros(rows, cols, 3, numFrames, 'uint8');
encodedTotalBytes = 0;

% Read and convert to YCbCr
originalYCbCrFrames = zeros(rows, cols, 3, numFrames, 'double');
for i = 1:numFrames
    frameRGB = read(vidObj, i);
    originalYCbCrFrames(:,:,:,i) = rgb2ycbcr(frameRGB);
end

% Main loop
for i = 1:numFrames
    rawFrame = originalYCbCrFrames(:,:,:,i);
    frameType = frameTypes{i};

    Y = rawFrame(:,:,1);
    Cb = rawFrame(:,:,2);
    Cr = rawFrame(:,:,3);

    if frameType == 'I'
        encodedY = jpeg_encode(Y, q_intra);
        encodedCb = jpeg_encode(Cb, q_intra);
        encodedCr = jpeg_encode(Cr, q_intra);

        decodedY = jpeg_decode(encodedY, q_intra);
        decodedCb = jpeg_decode(encodedCb, q_intra);
        decodedCr = jpeg_decode(encodedCr, q_intra);

        encodedFrames{i} = struct('Type', 'I', 'Y', encodedY, 'Cb', encodedCb, 'Cr', encodedCr);
        frameDecoded = cat(3, decodedY, decodedCb, decodedCr);
        decodedRGBFrames(:,:,:,i) = ycbcr2rgb(uint8(min(max(frameDecoded, 0), 255)));

        encodedTotalBytes = encodedTotalBytes + (nnz(encodedY) + nnz(encodedCb) + nnz(encodedCr)) * 2;

    elseif frameType == 'P'
        prevIdx = find(strcmp(frameTypes(1:i-1), 'I') | strcmp(frameTypes(1:i-1), 'P'), 1, 'last');
        if ~isempty(prevIdx)
            refY = double(decodedRGBFrames(:,:,1,prevIdx));
            refCb = double(decodedRGBFrames(:,:,2,prevIdx));
            refCr = double(decodedRGBFrames(:,:,3,prevIdx));

            diffY = motion_compensate(Y, refY, blockSize, searchRange, skipThreshold);
            diffCb = motion_compensate(Cb, refCb, blockSize, searchRange, skipThreshold);
            diffCr = motion_compensate(Cr, refCr, blockSize, searchRange, skipThreshold);

            encodedY = jpeg_encode(diffY, q_inter);
            encodedCb = jpeg_encode(diffCb, q_inter);
            encodedCr = jpeg_encode(diffCr, q_inter);

            decodedY = jpeg_decode(encodedY, q_inter) + refY;
            decodedCb = jpeg_decode(encodedCb, q_inter) + refCb;
            decodedCr = jpeg_decode(encodedCr, q_inter) + refCr;

            frameDecoded = cat(3, decodedY, decodedCb, decodedCr);
            decodedRGBFrames(:,:,:,i) = ycbcr2rgb(uint8(min(max(frameDecoded, 0), 255)));

            encodedFrames{i} = struct('Type', 'P', 'Y', encodedY, 'Cb', encodedCb, 'Cr', encodedCr, 'Ref', prevIdx);
            encodedTotalBytes = encodedTotalBytes + (nnz(encodedY) + nnz(encodedCb) + nnz(encodedCr)) * 2;
        end

    elseif frameType == 'B'
        prevIdx = find(strcmp(frameTypes(1:i-1), 'I') | strcmp(frameTypes(1:i-1), 'P'), 1, 'last');
        nextIdx = find(strcmp(frameTypes(i+1:end), 'I') | strcmp(frameTypes(i+1:end), 'P'), 1, 'first');
        if ~isempty(nextIdx)
            nextIdx = nextIdx + i;
        end
        if ~isempty(prevIdx) && ~isempty(nextIdx) && nextIdx <= numFrames
            refY = (double(decodedRGBFrames(:,:,1,prevIdx)) + double(decodedRGBFrames(:,:,1,nextIdx))) / 2;
            refCb = (double(decodedRGBFrames(:,:,2,prevIdx)) + double(decodedRGBFrames(:,:,2,nextIdx))) / 2;
            refCr = (double(decodedRGBFrames(:,:,3,prevIdx)) + double(decodedRGBFrames(:,:,3,nextIdx))) / 2;

            diffY = motion_compensate(Y, refY, blockSize, searchRange, skipThreshold);
            diffCb = motion_compensate(Cb, refCb, blockSize, searchRange, skipThreshold);
            diffCr = motion_compensate(Cr, refCr, blockSize, searchRange, skipThreshold);

            encodedY = jpeg_encode(diffY, q_inter);
            encodedCb = jpeg_encode(diffCb, q_inter);
            encodedCr = jpeg_encode(diffCr, q_inter);

            decodedY = jpeg_decode(encodedY, q_inter) + refY;
            decodedCb = jpeg_decode(encodedCb, q_inter) + refCb;
            decodedCr = jpeg_decode(encodedCr, q_inter) + refCr;

            frameDecoded = cat(3, decodedY, decodedCb, decodedCr);
            decodedRGBFrames(:,:,:,i) = ycbcr2rgb(uint8(min(max(frameDecoded, 0), 255)));

            encodedFrames{i} = struct('Type', 'B', 'Y', encodedY, 'Cb', encodedCb, 'Cr', encodedCr, 'RefPrev', prevIdx, 'RefNext', nextIdx);
            encodedTotalBytes = encodedTotalBytes + (nnz(encodedY) + nnz(encodedCb) + nnz(encodedCr)) * 2;
        end
    end
end

% Save original video (converted to RGB for visualization)
vOrig = VideoWriter('original_video.avi');
open(vOrig)
for i = 1:numFrames
    writeVideo(vOrig, ycbcr2rgb(uint8(originalYCbCrFrames(:,:,:,i))));
end
close(vOrig)

% Save decoded video
vDec = VideoWriter('decoded_video.avi');
open(vDec)
for i = 1:numFrames
    writeVideo(vDec, decodedRGBFrames(:,:,:,i));
end
close(vDec)

% Compute metrics
rawSize = rows * cols * 3 * numFrames;
rawSizeKB = rawSize / 1024;
encodedSizeKB = encodedTotalBytes / 1024;
compressionRatio = rawSizeKB / encodedSizeKB;
mse = mean((double(decodedRGBFrames(:)) - double(uint8(originalYCbCrFrames(:)))) .^ 2);
psnrVal = 10 * log10(255^2 / mse);

% Display results
fprintf('Raw Size: %.2f KB\n', rawSizeKB);
fprintf('Encoded Size: %.2f KB\n', encodedSizeKB);
fprintf('Compression Ratio: %.2f\n', compressionRatio);
fprintf('MSE: %.2f\n', mse);
fprintf('PSNR: %.2f dB\n', psnrVal);

%% --- Motion Compensation Function ---
function motionDiff = motion_compensate(target, ref, blockSize, range, threshold)
    [rows, cols] = size(target);
    motionDiff = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            i1 = min(i+blockSize-1, rows);
            j1 = min(j+blockSize-1, cols);
            block = double(target(i:i1, j:j1));
            bestMatch = ref(i:i1, j:j1);
            minError = inf;

            for dx = -range:range
                for dy = -range:range
                    x = i + dx; y = j + dy;
                    x1 = x + (i1 - i); y1 = y + (j1 - j);
                    if x >= 1 && y >= 1 && x1 <= rows && y1 <= cols
                        candidate = ref(x:x1, y:y1);
                        error = sum(abs(block(:) - candidate(:)));
                        if error < minError
                            minError = error;
                            bestMatch = candidate;
                        end
                    end
                end
            end
            if minError < threshold * numel(block)
                motionDiff(i:i1, j:j1) = 0;
            else
                motionDiff(i:i1, j:j1) = block - bestMatch;
            end
        end
    end
end

%% --- JPEG Encode ---
function encodedBlock = jpeg_encode(block, Q)
    blockSize = 8;
    [rows, cols] = size(block);
    encodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            sub = block(i:i+blockSize-1, j:j+blockSize-1);
            dctB = dct2(sub);
            quant = round(dctB ./ Q);
            encodedBlock(i:i+blockSize-1, j:j+blockSize-1) = quant;
        end
    end
end

%% --- JPEG Decode ---
function decodedBlock = jpeg_decode(encodedBlock, Q)
    blockSize = 8;
    [rows, cols] = size(encodedBlock);
    decodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            quant = encodedBlock(i:i+blockSize-1, j:j+blockSize-1);
            dctB = quant .* Q;
            idctB = idct2(dctB);
            decodedBlock(i:i+blockSize-1, j:j+blockSize-1) = idctB;
        end
    end
end
