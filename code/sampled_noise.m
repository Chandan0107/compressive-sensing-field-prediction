function [b] = sampled_noise(y,snr)
% Adds complex white Gaussian noise to each sample at the given SNR (dB),
% measured per sample. Same as awgn(y(i),snr,'measured') without needing
% the Communications Toolbox.
b = zeros(size(y));
for i = 1:length(y)
    noise_power = abs(y(i))^2 / 10^(snr/10);
    b(i) = y(i) + sqrt(noise_power/2) * (randn + 1j*randn);
end
end
