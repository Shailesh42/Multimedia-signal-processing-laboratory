% Video file path
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

% Frame types (I-frame, P-frame, B-frame)
frameTypes = repmat({'I', 'P', 'B', 'B', 'P', 'B', 'B'}, 1, ceil(numFrames / 7)); % Alternate frame types
frameTypes = frameTypes(1:numFrames); % Make sure we have the correct length for frame types

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

% Function for scalability techniques
function apply_scalability_technique(scalability_type, originalYCbCrFrames, q_intra, q_inter, numFrames, rows, cols, encodedFrames, decodedRGBFrames, frameTypes)
    encodedTotalBytes = 0;

    % Loop through all frames
    for i = 1:numFrames
        rawFrame = originalYCbCrFrames(:,:,:,i);
        frameType = frameTypes{i}; % Get the frame type ('I' or 'P')

        Y = rawFrame(:,:,1);   % Y (luminance) channel
        Cb = rawFrame(:,:,2);  % Cb (chrominance) channel
        Cr = rawFrame(:,:,3);  % Cr (chrominance) channel

        % Handle the spatial scalability technique
        if strcmp(scalability_type, 'spatial')
            % Spatial scalability: Downsample the frame by resizing
            if frameType == 'I'  % Intra-coded frame
                Y_resized = imresize(Y, 0.5);  % Downsample Y
                Cb_resized = imresize(Cb, 0.5); % Downsample Cb
                Cr_resized = imresize(Cr, 0.5); % Downsample Cr

                % Create the downsampled frame (Y, Cb, Cr)
                frameDecoded = cat(3, Y_resized, Cb_resized, Cr_resized);

                % Perform JPEG Encoding on the downsampled frame
                encodedY = jpeg_encode(Y_resized, q_intra);
                encodedCb = jpeg_encode(Cb_resized, q_intra);
                encodedCr = jpeg_encode(Cr_resized, q_intra);

                % Decode the frame back
                decodedY = jpeg_decode(encodedY, q_intra);
                decodedCb = jpeg_decode(encodedCb, q_intra);
                decodedCr = jpeg_decode(encodedCr, q_intra);

                % Reconstruct RGB frame from YCbCr
                decodedRGB = ycbcr2rgb(uint8(min(max(cat(3, decodedY, decodedCb, decodedCr), 0), 255)));

                % Assign the decoded RGB frame to the decodedRGBFrames array
                decodedRGBFrames(:,:,:,i) = decodedRGB; % This will now work because of consistent resizing

                % Store the encoded frames and compute the total size
                encodedFrames{i} = struct('Type', 'I', 'Y', encodedY, 'Cb', encodedCb, 'Cr', encodedCr);
                encodedTotalBytes = encodedTotalBytes + (numel(encodedY) + numel(encodedCb) + numel(encodedCr)) * 2;
            end

            % For predicted frames (P-frames)
        elseif frameType == 'P'
            % Find previous I or P frame for motion compensation
            prevIdx = find(strcmp(frameTypes(1:i-1), 'I') | strcmp(frameTypes(1:i-1), 'P'), 1, 'last');
            if ~isempty(prevIdx)
                refY = double(decodedRGBFrames(:,:,1,prevIdx)); % Reference frame Y
                refCb = double(decodedRGBFrames(:,:,2,prevIdx)); % Reference frame Cb
                refCr = double(decodedRGBFrames(:,:,3,prevIdx)); % Reference frame Cr

                % Perform motion compensation (this part needs to be defined in your motion_compensate function)
                diffY = motion_compensate(Y, refY); % Motion-compensated Y
                diffCb = motion_compensate(Cb, refCb); % Motion-compensated Cb
                diffCr = motion_compensate(Cr, refCr); % Motion-compensated Cr

                % Perform JPEG encoding for the motion-compensated frames
                encodedY = jpeg_encode(diffY, q_inter);
                encodedCb = jpeg_encode(diffCb, q_inter);
                encodedCr = jpeg_encode(diffCr, q_inter);

                % Decode the frames and add the reference frame (motion compensation)
                decodedY = jpeg_decode(encodedY, q_inter) + refY;
                decodedCb = jpeg_decode(encodedCb, q_inter) + refCb;
                decodedCr = jpeg_decode(encodedCr, q_inter) + refCr;

                % Reconstruct the full frame (Y, Cb, Cr)
                frameDecoded = cat(3, decodedY, decodedCb, decodedCr);

                % Convert to RGB
                decodedRGBFrames(:,:,:,i) = ycbcr2rgb(uint8(min(max(frameDecoded, 0), 255))); 

                % Store the encoded frames and compute the total size
                encodedFrames{i} = struct('Type', 'P', 'Y', encodedY, 'Cb', encodedCb, 'Cr', encodedCr, 'Ref', prevIdx);
                encodedTotalBytes = encodedTotalBytes + (numel(encodedY) + numel(encodedCb) + numel(encodedCr)) * 2;
            end
        end
    end

    % Save output video to file
    vDec = VideoWriter([scalability_type, '_decoded_video.avi']);
    open(vDec)
    for i = 1:numFrames
        writeVideo(vDec, decodedRGBFrames(:,:,:,i));
    end
    close(vDec)

    % Compute performance metrics
    mse = mean((double(decodedRGBFrames(:)) - double(uint8(originalYCbCrFrames(:)))) .^ 2);
    psnrVal = 10 * log10(255^2 / mse);
    ssimVal = ssim(uint8(decodedRGBFrames), uint8(originalYCbCrFrames));

    % Display results
    fprintf('%s Scalability:\n', scalability_type);
    fprintf('MSE: %.2f\n', mse);
    fprintf('PSNR: %.2f dB\n', psnrVal);
    fprintf('SSIM: %.2f\n', ssimVal);
    fprintf('Encoded Size: %.2f KB\n', encodedTotalBytes / 1024);
    fprintf('Compression Ratio: %.2f\n', (rows * cols * 3 * numFrames) / (encodedTotalBytes));
end
% Apply Scalability Techniques
apply_scalability_technique('spatial', originalYCbCrFrames, q_intra, q_inter, numFrames, rows, cols, encodedFrames, decodedRGBFrames, frameTypes);
apply_scalability_technique('temporal', originalYCbCrFrames, q_intra, q_inter, numFrames, rows, cols, encodedFrames, decodedRGBFrames, frameTypes);
apply_scalability_technique('snr', originalYCbCrFrames, q_intra, q_inter, numFrames, rows, cols, encodedFrames, decodedRGBFrames, frameTypes);

% --- Motion Compensation Function ---
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

% --- JPEG Encode ---
function encodedBlock = jpeg_encode(block, Q)
    blockSize = 8;
    [rows, cols] = size(block);
    encodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            % Ensure block does not exceed image boundaries
            i1 = min(i + blockSize - 1, rows);
            j1 = min(j + blockSize - 1, cols);
            sub = block(i:i1, j:j1);
            % JPEG DCT, Quantization and Zigzag
            encodedBlock(i:i1, j:j1) = dct2(sub) / Q;
        end
    end
end

% --- JPEG Decode ---
function decodedBlock = jpeg_decode(encodedBlock, Q)
    blockSize = 8;
    [rows, cols] = size(encodedBlock);
    decodedBlock = zeros(rows, cols);
    for i = 1:blockSize:rows
        for j = 1:blockSize:cols
            % Ensure block does not exceed image boundaries
            i1 = min(i + blockSize - 1, rows);
            j1 = min(j + blockSize - 1, cols);
            sub = encodedBlock(i:i1, j:j1);
            % JPEG Inverse DCT and Dequantization
            decodedBlock(i:i1, j:j1) = idct2(sub * Q);
        end
    end
end
