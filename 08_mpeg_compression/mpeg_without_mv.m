
videoFile = 'visiontraffic.avi';
vidObj = VideoReader(videoFile);
numFrames = floor(vidObj.Duration * vidObj.FrameRate);
rows = vidObj.Height;
cols = vidObj.Width;

originalFrames = zeros(rows, cols, numFrames, 'uint8');
for i = 1:numFrames
    originalFrames(:,:,i) = rgb2gray(read(vidObj, i));
end

q_intra = [8 16 19 22 26 27 29 34;
           16 16 22 24 27 29 34 37;
           19 22 26 27 29 34 34 38;
           22 22 26 27 29 34 37 40;
           22 26 27 29 32 35 40 48;
           26 27 29 32 35 40 48 58;
           26 27 29 34 38 46 56 89;
           27 29 35 38 46 56 69 83];
q_inter = 16 * ones(8, 8);
encodedFrames = cell(1, numFrames);
decodedFrames = zeros(rows, cols, numFrames, 'uint8');
frameTypes = repmat({'I','B','B','P','B','B','P','B','B'}, 1, ceil(numFrames/9));
frameTypes = frameTypes(1:numFrames);
encodedTotalBytes = 0;
for i = 1:numFrames
    rawFrame = originalFrames(:,:,i);
    frameType = frameTypes{i};
    if frameType == 'I'
        encoded = jpeg_encode(double(rawFrame), q_intra);
        encodedFrames{i} = struct('Type', 'I', 'Data', encoded);
        decodedFrames(:,:,i) = uint8(jpeg_decode(encoded, q_intra));
    elseif frameType == 'P'
        prevIdx = find(strcmp(frameTypes(1:i-1), 'I') | strcmp(frameTypes(1:i-1), 'P'), 1, 'last');
        if ~isempty(prevIdx)
            refFrame = decodedFrames(:,:,prevIdx);
            motionDiff = double(rawFrame) - double(refFrame);
            encoded = jpeg_encode(motionDiff, q_inter);
            encodedFrames{i} = struct('Type', 'P', 'Data', encoded, 'Ref', prevIdx);
            decodedFrames(:,:,i) = uint8(jpeg_decode(encoded, q_inter) + double(refFrame));
        end
    elseif frameType == 'B'
        prevIdx = find(strcmp(frameTypes(1:i-1), 'I') | strcmp(frameTypes(1:i-1), 'P'), 1, 'last');
        nextIdx = find(strcmp(frameTypes(i+1:end), 'I') | strcmp(frameTypes(i+1:end), 'P'), 1, 'first') + i;
        if ~isempty(prevIdx) && ~isempty(nextIdx) && nextIdx <= numFrames
            refPrev = decodedFrames(:,:,prevIdx);
            refNext = decodedFrames(:,:,nextIdx);
            avgRef = (double(refPrev) + double(refNext)) / 2;
            motionDiff = double(rawFrame) - avgRef;
            encoded = jpeg_encode(motionDiff, q_inter);
            encodedFrames{i} = struct('Type', 'B', 'Data', encoded, 'RefPrev', prevIdx, 'RefNext', nextIdx);
            decodedFrames(:,:,i) = uint8(jpeg_decode(encoded, q_inter) + avgRef);
        end
    end
    encodedTotalBytes = encodedTotalBytes + nnz(encoded) * 2;
end
vOrig = VideoWriter('original_video1.avi');
open(vOrig)
for i = 1:numFrames
    writeVideo(vOrig, repmat(originalFrames(:,:,i), [1 1 3]));
end
close(vOrig)

vDec = VideoWriter('decoded_video1.avi');
open(vDec)
for i = 1:numFrames
    writeVideo(vDec, repmat(decodedFrames(:,:,i), [1 1 3]));
end
close(vDec)
rawSize = rows * cols * numFrames; 
rawSizeKB = rawSize / 1024;
encodedSizeKB = encodedTotalBytes / 1024;
compressionRatio = rawSizeKB / encodedSizeKB;
mse = mean((double(originalFrames(:)) - double(decodedFrames(:))).^2);
fprintf('Original Raw Size: %.2f KB\n', rawSizeKB);
fprintf('Encoded Size: %.2f KB\n', encodedSizeKB);
fprintf('True Compression Ratio (Raw vs Encoded): %.2f\n', compressionRatio);
fprintf('MSE: %.2f\n', mse);

function encodedBlock = jpeg_encode(block, Q)
    blockSize = 8;
    [rows, cols] = size(block);
    encodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            subBlock = block(i:i+blockSize-1, j:j+blockSize-1);
            dctBlock = dct2(subBlock);
            quantBlock = round(dctBlock ./ Q);
            encodedBlock(i:i+blockSize-1, j:j+blockSize-1) = quantBlock;
        end
    end
end

function decodedBlock = jpeg_decode(encodedBlock, Q)
    blockSize = 8;
    [rows, cols] = size(encodedBlock);
    decodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            quantBlock = encodedBlock(i:i+blockSize-1, j:j+blockSize-1);
            dctBlock = quantBlock .* Q;
            idctBlock = idct2(dctBlock);
            decodedBlock(i:i+blockSize-1, j:j+blockSize-1) = idctBlock;
        end
    end
end
