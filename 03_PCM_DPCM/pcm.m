clc;
clear; 
close all;
load train.mat  
audiowrite('train.wav', y, Fs);  
mu = 255;  
bit_values = [4, 8, 12];
audio = zeros(length(y), 1); 
for i = 2:length(y)
    audio(i) = y(i) - y(i-1); 
end
figure;
subplot(4, 2, 1);
plot(audio(1:100));
title('Differentially Encoded Audio');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
compressed_audio = sign(audio) .* (log(1 + mu * abs(audio)) / log(1 + mu));
mse_compressed = mean((compressed_audio - audio).^2);
subplot(4, 2, 2);
plot(compressed_audio(1:100));
title(['Compressed Audio (\mu = ', num2str(mu),', MSE = ', num2str(mse_compressed), ')']);
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

for i = 1:3
    bitdepth = 2^(bit_values(i) - 1);
    quantized_audio = round(compressed_audio * (bitdepth - 1)) / (bitdepth - 1); 
    expanded_audio = sign(quantized_audio) .* ((1 + mu).^abs(quantized_audio) - 1) / mu;
    final = zeros(length(y), 1); 
    final(1) = y(1); 
    for j = 2:length(y)
        final(j) = expanded_audio(j) + final(j-1);
    end
 
    mse_expanded = mean((expanded_audio - audio).^2); 
    mse_final = mean((final - y).^2); 
    subplot(4, 2, i*2 + 1);
    plot(expanded_audio(1:100));
    title(['Expanded Audio (MSE = ', num2str(mse_expanded), ', Bits = ', num2str(bit_values(i)), ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    subplot(4, 2, i*2 + 2);
    plot(final(1:100));
    title(['Final Audio (MSE = ', num2str(mse_final), ', Bits = ', num2str(bit_values(i)), ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
end


load train.mat  
audiowrite('train.wav', y, Fs);  
mu = 255;  
audio = y ;
figure;
subplot(4, 2, 1);
plot(audio(1:100));
title('Original train.mat Audio');
xlabel('Sample Index');
ylabel('Amplitude');
bit_values = [ 4, 8, 12];
compressed_audio = sign(audio) .* (log(1 + mu * abs(audio)) / log(1 + mu));
mse_compressed = mean((compressed_audio - audio).^2);
subplot(4, 2, 2);
plot(compressed_audio(1:100));
title(['Compressed Audio (\mu = ', num2str(mu),', MSE = ', num2str(mse_compressed), ')']);
xlabel('Sample Index');
ylabel('Amplitude');
for i = 1: 3
    bitdepth = 2^(bit_values(i)-1);
    floor_r = floor(compressed_audio * bitdepth);
    quantized_audio = floor_r / bitdepth; 
    expanded_audio = sign(quantized_audio) .* ((1 + mu).^abs(quantized_audio) - 1) / mu;
    scaled_audio = expanded_audio * max(abs(audio));
    mse_expanded = mean((expanded_audio - audio).^2);
    mse_scaled = mean((scaled_audio - audio).^2);
    subplot(4, 2, i*2 +1);
    plot(expanded_audio(1:100));
    title(['Expanded Audio (MSE = ', num2str(mse_expanded),', Bits = ', num2str(bit_values(i)), ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
    subplot(4, 2, i*2 +2);
    plot(scaled_audio(1:100));
    title(['scaled Audio (MSE = ', num2str(mse_scaled),', Bits = ', num2str(bit_values(i)), ')']);
    xlabel('Sample Index');
    ylabel('Amplitude');
end
