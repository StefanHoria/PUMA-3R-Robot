%% DATE PUMA - aici modifici doar datele, apoi rulezi fișierul (Run / F5)
% Tot codul se află în PUMA_simulare.m (trebuie să fie în același folder).
% Unghiurile se dau în GRADE, lungimile în METRI.
% O etapă se dezactivează lăsând datele ei goale: []  (sau 0 la setări).
clear; clc; close all;

%% Dimensiunile robotului [m]
%            l3        l4       le
%         |------|---------|-------|
%  q1 ->  colonă l2, apoi cupla 2 (q2), cupla 3 (q3), efectorul E
date.l1 = 0.40;     % baza fixă, de la sol la cupla 1
date.l2 = 0.20;     % coloana verticală a elementului 1
date.l3 = 0.30;     % brațul orizontal al elementului 1 (până la cupla 2)
date.l4 = 0.40;     % elementul 2 (cupla 2 -> cupla 3)
date.le = 0.30;     % elementul 3 (cupla 3 -> efector)

%% Material și sarcină
date.raza = 0.03;        % raza elementelor [m]
date.rho  = 2700;        % densitatea [kg/m^3] (aluminiu)
date.mp   = 0.5;         % sarcina din efector [kg]
date.Fv   = [0 0 0];     % frecare vâscoasă în cuple [N·m·s/rad]

%% 1. Poziție statică [grade] - robotul desenat într-o singură configurație
%  q2 < 0 ridică brațul, q2 > 0 îl coboară (la fel q3)
date.q_static = [180  -90  90];          % [q1 q2 q3]

%% 2. Mișcare în spațiul cuplelor [grade] - puncte de trecere, câte o linie
%  Robotul trece pe rând prin fiecare configurație.
date.q_puncte = [   0   -30   60;
                   90   -60   90;
                  -45   -15   45;
                    0   -30   60 ];

%% 3. Urmărire r(t) [m] - puncte prin care trebuie să treacă efectorul
%  Între puncte efectorul merge în linie dreaptă. [] = etapa nu se rulează.
date.r_puncte = [ 0.75  -0.35  0.55;
                  0.50   0.60  0.85;
                  0.40   0.20  0.40 ];
date.cot = +1;           % +1 = cot sus, -1 = cot jos

%% Timp
date.T_segment = 3;      % durata unei mișcări între două puncte [s]
                         % (un număr, sau câte o valoare pentru fiecare segment)
date.N_segment = 100;    % numărul de puncte de calcul pe segment
date.grad      = 5;      % legea de mișcare: 3 = cubic, 5 = grad 5 (pornire lină)

%% 4. Mișcare liberă (tau = 0) - robotul cade sub greutatea proprie
date.miscare_libera = 1;     % 1 = da, 0 = nu  (pornește din q_static)
date.t_liber        = 2;     % durata [s]

%% 5. Acționări: motor de curent continuu + reductor planetar (câte unul pe cuplă)
%  Date de catalog maxon (pagina produsului, maxongroup.com); câte o valoare pentru cuplele 1, 2, 3.
%  Motor RE 40 Ø40 mm, 150 W, 24 V (148867); reductoare GP 52 C 230:1 (223101)
%  și GP 42 C 74:1 (203123) - www.maxongroup.com
date.motor.nume    = {'RE 40 148867 + GP 52 C 230:1 (223101)', ...
                      'RE 40 148867 + GP 52 C 230:1 (223101)', ...
                      'RE 40 148867 + GP 42 C 74:1 (203123)'};
date.motor.U_nom   = [24 24 24];              % tensiunea nominală [V]
date.motor.R       = [0.299 0.299 0.299];     % rezistența bobinajului [Ohm]
date.motor.km      = [30.2 30.2 30.2]*1e-3;   % constanta de cuplu [N·m/A]
date.motor.J_rotor = [142 142 142]*1e-7;      % inerția rotorului [kg·m^2] (142 g·cm^2)
date.motor.I0      = [137 137 137]*1e-3;      % curentul de mers în gol [A]
date.motor.n0      = [7580 7580 7580];        % turația de mers în gol [rpm]
date.motor.M_nom   = [177 177 177]*1e-3;      % cuplul nominal (max. continuu) [N·m]
date.motor.I_nom   = [6.0 6.0 6.0];           % curentul nominal (max. continuu) [A]
date.motor.ir      = [230 230 74];            % raportul de transmisie al reductorului
date.motor.ran     = [0.68 0.68 0.72];        % randamentul reductorului
date.motor.J_red   = [16.8 16.8 15]*1e-7;     % inerția reductorului (la intrare) [kg·m^2]
date.motor.M_red   = [30 30 15];              % cuplul max. continuu al reductorului [N·m]
date.motor.n_red   = [6000 6000 8000];        % turația max. de intrare a reductorului [rpm]

%% 6. Control PI (structura modelului prof. Pozna, realizată în Simulink)
date.control      = 1;          % 1 = da, 0 = nu
date.control_ref  = 'r';        % 'r' = urmărirea r(t) din r_puncte, 'q' = mișcarea prin q_puncte
date.poli         = [1/3 1/3];  % p1, p2 ca fracțiune din k1 (p3 = k1 - p1 - p2 rezultă)
date.acordare     = 'medie';    % inerția brațului la acordare: 'medie', 'max' sau 'motor' (neglijată)
date.t_mentinere  = 1;          % timp de menținere după terminarea mișcării [s]
date.simulink     = 1;          % 1 = construiește și rulează modelul Simulink contrPUMA.slx
date.rupere_bucla = 'filtru';   % 'filtru' (Tf mic) sau 'memory' (ca la Pozna, doar dacă k3*M_ii < 1)
date.T_filtru     = 1e-3;       % constanta de timp a filtrului de pe Tau [s]
date.deschide_simulink = 1;     % 1 = deschide modelul în Simulink după simulare

%% Ce se afișează
date.grafice      = 1;   % graficele de poziție/viteză/accelerație
date.dinamica     = 1;   % graficele cuplurilor
date.animatie     = 1;   % animația robotului
date.verificari   = 1;   % verificările modelului (în Command Window)
date.viteza_animatie = 1; % 1 = timp real, 2 = de 2 ori mai repede, 0.5 = cu încetinitorul

%% Rulare
rez = PUMA_simulare(date);
