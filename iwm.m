function denoised_y=iwm(time,input_singal,level)
% Decompose the signal using the wavelet transform

w = 'haar';
levels = level;
[C, L] = wavedec(input_singal, levels, w);

% Estimate the variance of each wavelet coefficient
va = zeros(1, levels+1);
for i = 1:levels+1
    va(i) = var(C((L(i)+1):(L(i+1))));
end

% Threshold the wavelet coefficients based on their variance
threshold = sqrt(2*log(length(input_singal))) * sqrt(va);
for i = 2:levels+1
    C((L(i-1)+1):(L(i))) = wthresh(C((L(i-1)+1):(L(i))), 's', threshold(i));
end

% Reconstruct the signal from the denoised wavelet coefficients
denoised_y = waverec(C, L, w);

% Plot the original and denoised signals

stairs(time, input_singal, 'b', 'LineWidth', 1)
hold on
stairs(time, denoised_y, 'r', 'LineWidth', 1)
legend('Original signal', 'Denoised signal')
xlabel('Time')
ylabel('Date Rate')
title('Independent Wavelet Model Denoising')
