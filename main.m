%% Donnees
rho = 988;              % kg/m3
mu  = 0.54e-3;          % Pa.s
L   = 40;               % m
eps = 0.1e-3;           % m
D   = 20e-3;            % m
Ks  = 6*0.5 + 2.5;      % singularites
DP_emetteur = 10000;    % Pa
eta_global = 0.65;

%% ETAPE 1 : courbe reseau
q_h = linspace(0, 2, 200);      % debit en m3/h
q_si = q_h / 3600;              % m3/s

A = pi*D^2/4;

v = zeros(size(q_h));
Re = zeros(size(q_h));
lambda = zeros(size(q_h));
DP_reseau = zeros(size(q_h));

for i = 1:length(q_h)
    if q_si(i) == 0
        v(i) = 0;
        Re(i) = 0;
        lambda(i) = 0;
        DP_reseau(i) = DP_emetteur;
    else
        v(i) = q_si(i)/A;
        Re(i) = rho*v(i)*D/mu;

        if Re(i) < 2300
            lambda(i) = 64/Re(i);
        else
            lambda(i) = 1 / (-1.8*log10(((eps/D)/3.7)^1.11 + 6.9/Re(i)))^2;
        end

        DP_L = lambda(i)*(L/D)*(rho*v(i)^2/2);
        DP_S = Ks*(rho*v(i)^2/2);

        DP_reseau(i) = DP_L + DP_S + DP_emetteur;
    end
end

%% ETAPE 2 : courbe pompe et point de fonctionnement
% IMPORTANT :
% ici q_h est pris en m3/h pour que la loi de pompe soit coherente
DP_pompe = -20000*q_h.^2 + 60000;

% Point de fonctionnement = minimum de l'ecart absolu
[~, idx] = min(abs(DP_reseau - DP_pompe));

q_star_h = q_h(idx);            % m3/h
q_star_si = q_star_h/3600;      % m3/s
DP_star = DP_reseau(idx);       % Pa
H_star = DP_star / (rho*9.81);  % mCE

fprintf('--- Point de fonctionnement ---\n');
fprintf('q* = %.3f m3/h\n', q_star_h);
fprintf('DP* = %.1f Pa\n', DP_star);
fprintf('H* = %.3f mCE\n', H_star);

figure;
plot(q_h, DP_reseau, 'b', 'LineWidth', 2); hold on;
plot(q_h, DP_pompe, 'r', 'LineWidth', 2);
plot(q_star_h, DP_star, 'ko', 'MarkerFaceColor', 'k');
grid on;
xlabel('Debit q_v (m^3/h)');
ylabel('\DeltaP (Pa)');
title('Courbes reseau / pompe');
legend('Reseau', 'Pompe', 'Point de fonctionnement', 'Location', 'best');

%% ETAPE 3 : comparaison aeraulique
Pth = 10000;           % W
cp_air = 1005;         % J/kg.K
DT_air = 15;           % degC
rho_air = 1.2;         % kg/m3 (hypothese usuelle)

mdot_air = Pth/(cp_air*DT_air);
qv_air = mdot_air/rho_air;      % m3/s

v_air = 5;                      % m/s
D_air = sqrt(4*qv_air/(pi*v_air));

fprintf('\n--- Comparaison air ---\n');
fprintf('Debit d air = %.4f m3/s = %.1f m3/h\n', qv_air, qv_air*3600);
fprintf('Diametre gaine air = %.3f m = %.1f mm\n', D_air, D_air*1000);

% Puissance hydraulique eau
P_hyd_eau = DP_star * q_star_si;
P_elec_eau = P_hyd_eau / eta_global;

fprintf('Puissance hydraulique eau = %.2f W\n', P_hyd_eau);
fprintf('Puissance electrique eau = %.2f W\n', P_elec_eau);

%% ETAPE 4 : cout annuel
% Remplacer selon ton hypothese
nb_heures = 2500;      % h/an
prix_kWh = 0.20;       % €/kWh

E_annuelle = (P_elec_eau/1000) * nb_heures;   % kWh/an
cout_annuel = E_annuelle * prix_kWh;

fprintf('\n--- Cout annuel ---\n');
fprintf('Energie annuelle = %.2f kWh/an\n', E_annuelle);
fprintf('Cout annuel = %.2f €/an\n', cout_annuel);

%% ETAPE 5 : optimisation des diametres
D_list_mm = [10 15 20 25 32];
q_test_h = 2.0;                  % debit de test en m3/h
q_test_si = q_test_h / 3600;
Pmax = 100;                      % W

Pelec_list = zeros(size(D_list_mm));
v_list = zeros(size(D_list_mm));
rejete = false(size(D_list_mm));

for i = 1:length(D_list_mm)
    D_i = D_list_mm(i)/1000;

    [DP_i, v_i] = calculPertes(q_test_si, D_i, L, eps, Ks, rho, mu, DP_emetteur);
    P_hyd_i = DP_i * q_test_si;
    P_elec_i = P_hyd_i / eta_global;

    Pelec_list(i) = P_elec_i;
    v_list(i) = v_i;

    if (v_i > 2) || (P_elec_i > Pmax)
        rejete(i) = true;
    end
end

fprintf('\n--- Optimisation diametres ---\n');
for i = 1:length(D_list_mm)
    etat = "ACCEPTE";
    if rejete(i)
        etat = "REJETE";
    end
    fprintf('D = %2d mm | v = %.3f m/s | P_elec = %.2f W | %s\n', ...
        D_list_mm(i), v_list(i), Pelec_list(i), etat);
end

figure;
bar(D_list_mm, Pelec_list); hold on;
yline(Pmax, 'r-', 'LineWidth', 2);

idx_rej = find(rejete);
plot(D_list_mm(idx_rej), Pelec_list(idx_rej), 'ro', ...
    'MarkerSize', 10, 'LineWidth', 2);

grid on;
xlabel('Diametre (mm)');
ylabel('Puissance electrique (W)');
title('Optimisation : P_{elec} en fonction du diametre');
legend('P_{elec}', 'Seuil 100 W', 'Diametres rejetes', 'Location', 'best');
