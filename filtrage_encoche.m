function y = filtrage_encoche(x,f0,fe,epsilon)
% =========================================================================
% Projet ECG 2014
% Quentin Biache
% Anthony Delannoy
%
% La fonction filtrage_encoche réalise le filtrage suivant la formule
% donnée dans le sujet de TP avec les paramètres :
%   - x : signal à filtrer
%   - f0 : la fréquence à filtrer
%   - fe : la fréquence d'échantillonnage
%   - epsilon : la sélectivité du filtre
% =========================================================================

w0 = 2*pi*f0/fe;
a1 = -2*cos(w0);

% Calcul des coefficients
B = [1 a1 1];
A = [1 (1-epsilon)*a1 (1-epsilon).^2];

y = filter(B,A,x);
end

