function y = filtrage_encoche(x,f0,fe,epsilon)
% =========================================================================
% Projet ECG 2014
% Quentin Biache
% Anthony Delannoy
%
% La fonction filtrage_encoche r�alise le filtrage suivant la formule
% donn�e dans le sujet de TP avec les param�tres :
%   - x : signal � filtrer
%   - f0 : la fr�quence � filtrer
%   - fe : la fr�quence d'�chantillonnage
%   - epsilon : la s�lectivit� du filtre
% =========================================================================

w0 = 2*pi*f0/fe;
a1 = -2*cos(w0);

% Calcul des coefficients
B = [1 a1 1];
A = [1 (1-epsilon)*a1 (1-epsilon).^2];

y = filter(B,A,x);
end

