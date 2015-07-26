function [autocorr_biased, autocorr_unbiased, x_corr_biased, x_corr_unbiased] = calcul_autocorr(s, N, N_fft)
% =========================================================================
% Projet ECG 2014
% Quentin Biache
% Anthony Delannoy
%
% La fonction calcul_autocorr se charge de renvoyer tous les paramètres
% nécessaires pour calculer un corrélogramme biaisé ou non biaisé avec zero
% padding avec les arguments :
%   - s : le signal 
%   - N : le nombre de points de signal
%   - N_fft : le nombre de points avec zero padding
% =========================================================================



autocorr_biased = xcorr(s,'biased');
% Ajout de zeros pour le zero padding
autocorr_biased = [autocorr_biased(N:2*N-1,1);zeros(N_fft-length(autocorr_biased),1);autocorr_biased(2*N-1:-1:N+1,1)];
x_corr_biased = abs(fft(autocorr_biased));

autocorr_unbiased = xcorr(s,'unbiased');
autocorr_unbiased = [autocorr_unbiased(N:2*N-1,1);zeros(N_fft-length(autocorr_unbiased),1);autocorr_unbiased(2*N-1:-1:N+1,1)];
x_corr_unbiased = abs(fft(autocorr_unbiased));
end

