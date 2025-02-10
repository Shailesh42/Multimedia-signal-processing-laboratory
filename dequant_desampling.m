clc;
close all;
load handel.mat
audiowrite('handel.wav', y, Fs);
sound(y, Fs);

var = [0.02, 0.05];
filtersize = [3, 5, 11, 13];                            
wightfilterkernal{1} = [1, 2, 1];                     
wightfilterkernal{2} = [1, 2, 4,2, 1];                 
wightfilterkernal{3} = [1, 2, 4, 8, 16, 32, 16, 8, 4, 2, 1];   
wightfilterkernal{4} = [1, 2, 4, 8, 16, 32, 64, 32, 16, 8, 4, 2, 1]; 
n = length(y);
a = audioread("handel.wav");
figure();
subplot(3,1,1);
plot(a(1:500));
title('original Audio');
xlabel('Time (s)');     ylabel('Amplitude');
noise = sqrt(var(1)) * randn(size(y));
subplot(3, 1, 2);
audio = y+ noise ; plot(audio(1:500));
title(['Noisy Audio (Variance = ', num2str(var(1)), ')']);
xlabel('Time (s)');     ylabel('Amplitude');
noise = sqrt(var(2)) * randn(size(y));
subplot(3, 1, 3);
audio = y+ noise ; plot(audio(1:500));
title(['Noisy Audio (Variance = ', num2str(var(2)), ')']);
xlabel('Time (s)');     ylabel('Amplitude');
for i = 1:2
    noise = sqrt(var(i)) * randn(size(y));
    audio = y + noise;
    mse_values = zeros(1, 4);
    figure();
    for j = 1:4
        start = filtersize(j);                
        hole = floor(start / 2);               
        y_mean = zeros(size(y));               
        for k = hole + 1 : n - hole
            window = audio(k - hole : k + hole);
            y_mean(k) = mean(window);
        end
        mse_values(j) = mean((y - y_mean).^2);             
        subplot(4, 1, j );
        plot(y_mean(1:500));
        title(['Mean Filter (Variance = ', num2str(var(i)), ', FilterSize = ', num2str(filtersize(j)), ', MSE = ', num2str(mse_values(j), '%.4f'), ')']);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    disp(['Variance = ', num2str(var(i)), ' (Mean Filter Results)']);
    disp('MSE values:');
    disp(mse_values);
   
    figure();
    for j = 1:4
        start = filtersize(j);                
        hole = floor(start / 2);               
        y_median = zeros(size(y));             
        for k = hole + 1 : n - hole
            window = audio(k - hole : k + hole);
            y_median(k) = median(window);
        end
        mse_values(j) = mean((y - y_median).^2);
        subplot(4, 1, j );
        plot(y_median(1:500));
        title(['Median Filter (Variance = ', num2str(var(i)), ', FilterSize = ', num2str(filtersize(j)), ', MSE = ', num2str(mse_values(j), '%.4f'), ')']);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    disp(['Variance = ', num2str(var(i)), ' (Median Filter Results)']);
    disp('MSE values:');
    disp(mse_values);
    
    figure();
    for j = 1:4
        kernel = wightfilterkernal{j};             % Select kernel for current filter
        kernel_len = length(kernel);               % Length of the kernel
        half_len = floor(kernel_len / 2);
        y_weight = zeros(size(y)); 
        for k = half_len + 1 : n - half_len
            window = audio(k - half_len : k + half_len);  
            weighted_sum = 0;
            for m = 1:kernel_len
                 weighted_sum = weighted_sum + window(m) * kernel(m);  
            end
            y_weight(k) = weighted_sum / sum(kernel);  
        end
        mse_values(j) = mean((y - y_weight).^2);            
        subplot(4, 1, j);
        plot(y_weight(1:500));
        title(['Weight Filter (Variance = ', num2str(var(i)), ', FilterSize = ', num2str(filtersize(j)), ', MSE = ', num2str(mse_values(j), '%.4f'), ')']);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end  
    disp(['Variance = ', num2str(var(i)), ' (Weighted Filter Results)']);
    disp('MSE values:');
    disp(mse_values);
end
