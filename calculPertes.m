%% Fonction de calcul des pertes de charge
function [DP_tot, v] = calculPertes(qv, D, L, eps, Ks, rho, mu, DP_emetteur)
% qv en m3/s
% D en m
% retourne la perte de charge totale en Pa et la vitesse en m/s

A = pi*D^2/4;
v = qv / A;

if qv == 0
    DP_tot = DP_emetteur;
    return;
end

Re = rho*v*D/mu;

if Re < 2300
    lambda = 64/Re;
else
    lambda = 1 / (-1.8*log10(((eps/D)/3.7)^1.11 + 6.9/Re))^2;
end

DP_L = lambda*(L/D)*(rho*v^2/2);
DP_S = Ks*(rho*v^2/2);

DP_tot = DP_L + DP_S + DP_emetteur;
end
