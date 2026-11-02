function rez = PUMA_simulare(date)
% PUMA_SIMULARE  Simularea completă a structurii de poziționare PUMA (3 cuple R)
%
%   rez = PUMA_simulare(date)
%
%   Datele (dimensiuni, unghiuri, puncte, setări) se scriu în Date_PUMA.m.
%   Orice câmp care lipsește primește o valoare implicită.
%   Această funcție conține tot codul: modelul geometric (direct și invers),
%   modelul cinematic, modelul dinamic, vizualizarea și animația.
%
%   Etape (fiecare se poate opri lăsând datele ei goale [] sau 0):
%     1. poziția statică q_static        -> desen + 0_eT + cupluri gravitaționale
%     2. mișcare prin punctele q_puncte  -> animație, cinematică, cupluri
%     3. urmărire r(t) prin r_puncte     -> cinematică inversă, animație, cupluri
%     4. mișcare liberă (tau = 0)        -> dinamică directă, energie, animație
%
%   rez - structură cu toate rezultatele (t, Q, dQ, ddQ, P, V, A, Tau, ...)

if nargin < 1, date = struct(); end
date = completeaza_date(date);
par  = parametri_PUMA(date);
rez.par = par;

fprintf('=== Structura PUMA: l1=%.3f  l2=%.3f  l3=%.3f  l4=%.3f  le=%.3f m ===\n', ...
        par.l1, par.l2, par.l3, par.l4, par.le);
fprintf('Mase: element 2 = %.3f kg, element 3 + sarcină = %.3f kg\n', par.m2, par.m3);

%% 1. Poziția statică
if ~isempty(date.q_static)
    q = deg2rad(date.q_static(:));
    [p, T0e] = geometric(q, par);
    [M, ~, G] = model_dinamic(q, [], par);
    fprintf('\n--- 1. Poziție statică q = [%.1f  %.1f  %.1f] grade ---\n', date.q_static);
    disp('0_eT = '); disp(T0e);
    fprintf('Efector E: x = %.4f  y = %.4f  z = %.4f m\n', p);
    fprintf('Cupluri necesare pentru menținere (gravitație): tau = [%.3f  %.3f  %.3f] N·m\n', G);
    if date.verificari
        try
            q_inv = geometric_invers(p, par, cot_din_q(q));
            fprintf('Verificare MGD -> MGI -> MGD: eroare %.2e m\n', norm(geometric(q_inv, par) - p));
        catch
            fprintf('Verificare MGI: configurație în care efectorul trece de axa z0 (nu se verifică).\n');
        end
    end
    vizualizare_robot(q, par, [], '1. Poziție statică');
    rez.static = struct('q', q, 'p', p, 'T0e', T0e, 'M', M, 'G', G);
end

%% 2. Mișcare în spațiul cuplelor (prin punctele q_puncte)
if size(date.q_puncte, 1) >= 2
    fprintf('\n--- 2. Mișcare în spațiul cuplelor (%d puncte) ---\n', size(date.q_puncte, 1));
    [t, pvc] = traiectorie_puncte(deg2rad(date.q_puncte).', date.T_segment, date.N_segment, date.grad);
    Q = pvc(1:3,:);  dQ = pvc(4:6,:);  ddQ = pvc(7:9,:);

    et = '';  if date.grafice, et = '2. Spațiul cuplelor'; end
    grafic_cuple(t, Q, dQ, ddQ, et, 'Mișcare în spațiul cuplelor');
    [P, V, A, W] = cinematica_directa(Q, dQ, ddQ, t, par, et);
    if ~date.dinamica, et = ''; end
    Tau = dinamica_inversa(Q, dQ, ddQ, t, par, et);
    fprintf('Viteza maximă a efectorului: %.3f m/s\n', max(sqrt(sum(V.^2, 1))));
    fprintf('Cupluri maxime: |tau| = [%.3f  %.3f  %.3f] N·m\n', max(abs(Tau), [], 2));
    if date.verificari
        V_num = gradient(P, t(2) - t(1));
        fprintf('Verificare Jacobian (J*dq vs. dp/dt numeric): %.2e m/s\n', ...
                max(max(abs(V(:,2:end-1) - V_num(:,2:end-1)))));
    end
    if date.animatie
        animatie(Q, par, [], date.viteza_animatie, t, '2. Spațiul cuplelor');
    end
    rez.cuple = struct('t', t, 'Q', Q, 'dQ', dQ, 'ddQ', ddQ, 'P', P, 'V', V, 'A', A, 'W', W, 'Tau', Tau);
end

%% 3. Urmărirea traiectoriei r(t) (prin punctele r_puncte)
if size(date.r_puncte, 1) >= 2
    fprintf('\n--- 3. Urmărire r(t) (%d puncte) ---\n', size(date.r_puncte, 1));
    [t, rvc] = traiectorie_puncte(date.r_puncte.', date.T_segment, date.N_segment, date.grad);
    R = rvc(1:3,:);  dR = rvc(4:6,:);  ddR = rvc(7:9,:);
    verifica_spatiu_lucru(R, t, par);

    [Qd, dQd, ddQd] = cinematica_inversa(R, dR, ddR, par, date.cot);
    et = '';  if date.grafice, et = '3. Urmărire r(t)'; end
    grafic_cuple(t, Qd, dQd, ddQd, et, 'Cuple pentru urmărirea r(t)');
    [Pc, Vc, Ac] = cinematica_directa(Qd, dQd, ddQd, t, par, et);
    if ~date.dinamica, et = ''; end
    Tau = dinamica_inversa(Qd, dQd, ddQd, t, par, et);
    fprintf('Cupluri maxime: |tau| = [%.3f  %.3f  %.3f] N·m\n', max(abs(Tau), [], 2));
    if date.verificari
        fprintf('Verificare r(t): eroare poziție %.2e m, viteză %.2e m/s, accelerație %.2e m/s^2\n', ...
                max(abs(Pc(:) - R(:))), max(abs(Vc(:) - dR(:))), max(abs(Ac(:) - ddR(:))));
        e_dd = 0;
        for i = 1:numel(t)
            [Mi, Ci, Gi] = model_dinamic(Qd(:,i), dQd(:,i), par);
            e_dd = max(e_dd, norm(Mi \ (Tau(:,i) - Ci*dQd(:,i) - Gi - par.Fv(:).*dQd(:,i)) - ddQd(:,i)));
        end
        fprintf('Verificare dinamică directă(dinamică inversă): eroare %.2e rad/s^2\n', e_dd);
    end
    if date.animatie
        animatie(Qd, par, R, date.viteza_animatie, t, '3. Urmărire r(t)');
    end
    rez.urmarire = struct('t', t, 'R', R, 'dR', dR, 'ddR', ddR, 'Q', Qd, 'dQ', dQd, 'ddQ', ddQd, 'Tau', Tau);
end

%% 4. Mișcare liberă (tau = 0): robotul cade sub greutatea proprie
if date.miscare_libera
    if ~isempty(date.q_static)
        q0 = deg2rad(date.q_static(:));
    elseif ~isempty(date.q_puncte)
        q0 = deg2rad(date.q_puncte(1,:).');
    else
        q0 = [0; 0; 0];
    end
    fprintf('\n--- 4. Mișcare liberă din q = [%.1f  %.1f  %.1f] grade, %.1f s ---\n', rad2deg(q0), date.t_liber);
    t = linspace(0, date.t_liber, max(2, round(100*date.t_liber)) + 1);
    [Ql, dQl] = simulare_dinamica([q0; 0; 0; 0], t, @(tt, x) [0; 0; 0], par);
    E = zeros(size(t));
    for i = 1:numel(t)
        [Mi, ~, ~, Vi] = model_dinamic(Ql(:,i), dQl(:,i), par);
        E(i) = 0.5*dQl(:,i).'*Mi*dQl(:,i) + Vi;
    end
    if date.verificari && all(par.Fv == 0)
        fprintf('Verificare conservarea energiei: variație %.2e J (din %.3f J)\n', max(E) - min(E), E(1));
    end
    if date.grafice
        tl = tab_grafice('4. Mișcare liberă', 'Cuple și energie', 2, 1);
        ax = nexttile(tl); plot(ax, t, rad2deg(Ql), 'LineWidth', 1.5); grid(ax, 'on');
        ylabel(ax, 'q [°]'); legenda_tl(ax, {'q_1','q_2','q_3'}); title(ax, 'Mișcare liberă sub gravitație (\tau = 0)');
        ax = nexttile(tl); plot(ax, t, E - E(1), 'LineWidth', 1.5); grid(ax, 'on');
        ylabel(ax, 'E - E_0 [J]'); xlabel(ax, 't [s]'); title(ax, 'Variația energiei totale');
    end
    if date.animatie
        animatie(Ql, par, [], date.viteza_animatie, t, '4. Mișcare liberă');
    end
    rez.liber = struct('t', t, 'Q', Ql, 'dQ', dQl, 'E', E);
end

%% 5. Controlul PI al celor 3 motoare (structura din modelul prof. Pozna)
%  q_d(t) -> Controler (e = q_d - q -> PI -> saturație) -> v -> MotoReductor
%  (ddq = k2*v - k1*dq - k3*tau) -> q, dq, ddq -> StructuraPUMA (dinamica
%  inversă) -> tau -> înapoi la motoare (cuplul rezistent).
if date.control
    fprintf('\n--- 5. Control PI: motor cc + reductor pe fiecare cuplă ---\n');
    folder = fileparts(mfilename('fullpath'));

    % 5a. Traiectoria dorită a cuplelor (din r(t) prin modelul cinematic invers)
    [tc, Qc, Rc] = referinta_control(date, par);

    % 5b. Acționările și acordarea regulatoarelor PI prin alocarea polilor
    act = actionari(date, par, Qc);
    afiseaza_actionari(act);

    % 5c. Poz_dorita.mat - intrarea blocului From File (rândul 1 = timpul)
    Poz_dorita = [tc; Qc];
    save(fullfile(folder, 'Poz_dorita.mat'), 'Poz_dorita');
    t_stop = tc(end) + date.t_mentinere;

    % 5d. Simularea în MATLAB (aceeași buclă, cuplajul motor-structură rezolvat exact)
    res_m = control_matlab(tc, Qc, act, par, t_stop);
    res = res_m;  sursa = 'MATLAB (ode45)';

    % 5e. Modelul Simulink contrPUMA.slx (construit automat) și simularea lui
    if date.simulink
        if ~exist('new_system')                                    %#ok<EXIST>
            fprintf('Simulink nu este instalat - rămân rezultatele din MATLAB.\n');
        else
            fisier = fullfile(folder, [date.model_simulink '.slx']);
            if strcmpi(date.rupere_bucla, 'memory') && any(act.raport_inertii(:,3) >= 1)
                warning(['Blocul Memory (ca la Pozna) rupe bucla algebrică doar dacă k3*M_ii < 1; ' ...
                         'aici k3*M_ii ajunge la [%s] - simularea va diverge. Folosiți rupere_bucla = ''filtru''.'], ...
                        num2str(act.raport_inertii(:,3).', '%.2f  '));
            end
            res = control_simulink(fisier, act, par, Qc(:,1), t_stop, date, folder);
            sursa = sprintf('Simulink (%s.slx)', date.model_simulink);
            fprintf('Modelul Simulink a fost salvat: %s\n', fisier);
            if date.verificari
                Qm = interp1(res_m.t, res_m.Q.', res.t).';
                fprintf('Verificare Simulink vs. MATLAB: diferența maximă în cuple %.2e rad\n', ...
                        max(abs(Qm(:) - res.Q(:))));
            end
        end
    end

    % 5f. Verificarea acționărilor pe mișcarea simulată
    res = marimi_motor(res, tc, Qc, act);
    verifica_actionari(res, act);
    if ~isempty(Rc)
        Rd = interp1(tc, Rc.', min(res.t, tc(end))).';
        Pr = zeros(size(Rd));
        for i = 1:numel(res.t), Pr(:,i) = geometric(res.Q(:,i), par); end
        res.eroare_r = sqrt(sum((Pr - Rd).^2, 1));
        fprintf('Eroarea de urmărire a lui r(t) (efector): max %.3f mm\n', 1000*max(res.eroare_r));
    end
    if date.grafice
        grafice_control(res, tc, Qc, act, sursa);
    end
    if date.animatie
        tu = 0:0.02:res.t(end);
        animatie(interp1(res.t, res.Q.', tu).', par, Rc, date.viteza_animatie, tu, '5. Control PI');
    end
    rez.control = struct('act', act, 'tc', tc, 'Qd', Qc, 'Rd', Rc, 'res', res, ...
                         'res_matlab', res_m, 'sursa', sursa);
end
fprintf('\nGata.\n');
end


%% ========================================================================
%  FUNCȚII LOCALE
%  ========================================================================

%% Date implicite
function date = completeaza_date(date)
implicit = struct( ...
    'l1', 0.40, 'l2', 0.20, 'l3', 0.30, 'l4', 0.40, 'le', 0.30, ...
    'raza', 0.03, 'rho', 2700, 'mp', 0.5, 'g', 9.81, 'Fv', [0 0 0], ...
    'q_static', [180 -90 90], ...
    'q_puncte', [0 -30 60; 90 -60 90], ...
    'r_puncte', [], ...
    'cot', 1, 'T_segment', 3, 'N_segment', 100, 'grad', 5, ...
    'miscare_libera', 0, 't_liber', 2, ...
    'grafice', 1, 'dinamica', 1, 'animatie', 1, 'verificari', 1, 'viteza_animatie', 1, ...
    'control', 0, 'control_ref', 'r', 'dt_control', 1e-3, 't_mentinere', 1, ...
    'poli', [1/3 1/3], 'acordare', 'medie', 'simulink', 1, ...
    'rupere_bucla', 'filtru', 'T_filtru', 1e-3, 'deschide_simulink', 1, ...
    'model_simulink', 'contrPUMA');
campuri = fieldnames(implicit);
for k = 1:numel(campuri)
    if ~isfield(date, campuri{k})
        date.(campuri{k}) = implicit.(campuri{k});
    end
end
% Acționările implicite: maxon RE 40 (24 V) + reductoare planetare GP 52 C / GP 42 C
motor = struct( ...
    'nume',    {{'RE 40 148867 + GP 52 C 230:1 (223101)', ...
                 'RE 40 148867 + GP 52 C 230:1 (223101)', ...
                 'RE 40 148867 + GP 42 C 74:1 (203123)'}}, ...
    'U_nom',   [24 24 24], 'R', [0.299 0.299 0.299], 'km', [30.2 30.2 30.2]*1e-3, ...
    'J_rotor', [142 142 142]*1e-7, 'I0', [137 137 137]*1e-3, 'n0', [7580 7580 7580], ...
    'M_nom',   [177 177 177]*1e-3, 'I_nom', [6.0 6.0 6.0], ...
    'ir',      [230 230 74], 'ran', [0.68 0.68 0.72], 'J_red', [16.8 16.8 15]*1e-7, ...
    'M_red',   [30 30 15], 'n_red', [6000 6000 8000]);
if ~isfield(date, 'motor'), date.motor = struct(); end
campuri = fieldnames(motor);
for k = 1:numel(campuri)
    if ~isfield(date.motor, campuri{k})
        date.motor.(campuri{k}) = motor.(campuri{k});
    end
end
if numel(date.poli) ~= 2 || sum(date.poli) >= 1 || any(date.poli <= 0)
    error('poli = [p1 p2] (fracțiuni din k1) trebuie să fie pozitivi, cu p1 + p2 < 1.');
end
if ~isempty(date.q_puncte) && size(date.q_puncte, 2) ~= 3
    error('q_puncte trebuie să aibă 3 coloane: [q1 q2 q3] pe fiecare linie.');
end
if ~isempty(date.r_puncte) && size(date.r_puncte, 2) ~= 3
    error('r_puncte trebuie să aibă 3 coloane: [x y z] pe fiecare linie.');
end
end

function c = cot_din_q(q)
if sin(q(3)) < 0, c = -1; else, c = 1; end
end

%% Traiectorie prin mai multe puncte
function [t, pvc] = traiectorie_puncte(Pts, Tseg, Nseg, grad, dt)
% Pts - n x k: k puncte de trecere (coloane); între două puncte consecutive
% se generează o mișcare polinomială (viteză și accelerație nule în puncte).
% Tseg - durata fiecărui segment (scalar sau vector cu k-1 valori)
% dt   - (opțional) pas de timp fix; înlocuiește numărul de puncte Nseg
k = size(Pts, 2);
if isscalar(Tseg), Tseg = Tseg*ones(1, k-1); end
if numel(Tseg) ~= k-1
    error('T_segment trebuie să fie un număr sau un vector cu %d valori (câte una pe segment).', k-1);
end
t = []; pvc = []; t0 = 0;
for s = 1:k-1
    if nargin >= 5 && ~isempty(dt), Nseg = round(Tseg(s)/dt) + 1; end
    [ps, ts] = traiect(Pts(:, s:s+1), Tseg(s), Nseg, grad, 0);
    if s > 1, ps = ps(:, 2:end); ts = ts(2:end); end   % fără punct dublat
    t = [t, ts + t0];       %#ok<AGROW>
    pvc = [pvc, ps];        %#ok<AGROW>
    t0 = t0 + Tseg(s);
end
end

function verifica_spatiu_lucru(R, t, par)
% Verifică dacă toate punctele traiectoriei r(t) pot fi atinse
dmin = abs(par.l4 - par.le);  dmax = par.l4 + par.le;
for i = 1:size(R, 2)
    rho = sqrt(R(1,i)^2 + R(2,i)^2) - par.l3;
    h = par.l1 + par.l2 - R(3,i);
    d = sqrt(rho^2 + h^2);
    if d > dmax || d < dmin
        error(['Traiectoria r(t) iese din spațiul de lucru la t = %.2f s, punctul [%.3f %.3f %.3f] m.\n' ...
               'Distanța de la cupla 2 la efector este %.3f m, trebuie să fie între %.3f și %.3f m.\n' ...
               'Cupla 2 se află la înălțimea l1+l2 = %.3f m și la distanța l3 = %.3f m de axa z.'], ...
               t(i), R(:,i), d, dmin, dmax, par.l1 + par.l2, par.l3);
    end
end
end

function grafic_cuple(t, Q, dQ, ddQ, etapa, titlu)
% Poziția, viteza și accelerația cuplelor (etapa = '' -> fără grafic)
if isempty(etapa), return; end
tl = tab_grafice(etapa, 'Cuple (q, dq, ddq)', 3, 1);
g = {rad2deg(Q), 'q [°]'; rad2deg(dQ), 'dq [°/s]'; rad2deg(ddQ), 'ddq [°/s^2]'};
for k = 1:3
    ax = nexttile(tl);
    plot(ax, t, g{k,1}, 'LineWidth', 1.5); grid(ax, 'on'); ylabel(ax, g{k,2});
end
title(tl, titlu); xlabel(ax, 't [s]');
legenda_tl(ax, {'q_1', 'q_2', 'q_3'});
end

%% Parametrii robotului
function par = parametri_PUMA(date)
% Parametrii geometrici și inerțiali (din datele introduse)
%   {0} --[Rz(q1), t=(0,0,l1)]--> {1} --[Ry(q2), t=(l3,0,l2)]--> {2}
%       --[Ry(q3), t=(l4,0,0)]--> {3} --[I3,    t=(le,0,0)]--> {e}
% Elementele sunt cilindri omogeni de rază "raza" și densitate "rho".
par.l1 = date.l1;  par.l2 = date.l2;  par.l3 = date.l3;  par.l4 = date.l4;  par.le = date.le;
par.raza = date.raza;  par.rho = date.rho;  par.mp = date.mp;  par.g = date.g;
par.Fv = date.Fv(:);

r = par.raza;
masa = @(L) par.rho*pi*r^2*L;

% Element 1: doar inerția față de axa z0 contează
par.m1a = masa(par.l2);
par.m1b = masa(par.l3);
par.I1zz = par.m1a*r^2/2 + par.m1b*(3*r^2 + par.l3^2)/12 + par.m1b*(par.l3/2)^2;

% Element 2: bară l4 pe axa x2, CM la l4/2
par.m2  = masa(par.l4);
par.lc2 = par.l4/2;
par.I2  = [par.m2*r^2/2, par.m2*(3*r^2 + par.l4^2)/12, par.m2*(3*r^2 + par.l4^2)/12];  % [Ixx Iyy Izz] față de CM

% Element 3: bară le + sarcina mp în E
mb = masa(par.le);
par.m3  = mb + par.mp;
par.lc3 = (mb*par.le/2 + par.mp*par.le)/par.m3;
Iyy3 = mb*(3*r^2 + par.le^2)/12 + mb*(par.le/2 - par.lc3)^2 + par.mp*(par.le - par.lc3)^2;
par.I3  = [mb*r^2/2, Iyy3, Iyy3];
par.mb3 = mb;
end

%% Model geometric direct
function [p, T0e, T] = geometric(q, par, f)
% [p, T0e, T] = geometric(q, par, f)
%   f = 0 doar calcul, f = 1 deschide vizualizarea robotului
% Matricile din calcul:
%   0_1T = [Rz(q1) (0,0,l1)]   1_2T = [Ry(q2) (l3,0,l2)]
%   2_3T = [Ry(q3) (l4,0,0)]   3_eT = [I3     (le,0,0)]
%   0_2T = 0_1T*1_2T,  0_3T = 0_2T*2_3T,  0_eT = 0_3T*3_eT
if nargin < 3, f = 0; end
Rz = @(a) [cos(a) -sin(a) 0;  sin(a) cos(a) 0;  0 0 1];
Ry = @(a) [cos(a) 0 sin(a);  0 1 0;  -sin(a) 0 cos(a)];
omog = @(Rm, tv) [Rm, tv(:); 0 0 0 1];

T01 = omog(Rz(q(1)), [0; 0; par.l1]);
T12 = omog(Ry(q(2)), [par.l3; 0; par.l2]);
T23 = omog(Ry(q(3)), [par.l4; 0; 0]);
T3e = omog(eye(3),   [par.le; 0; 0]);

T02 = T01*T12;
T03 = T02*T23;
T0e = T03*T3e;

p = T0e(1:3, 4);
T = struct('T01', T01, 'T02', T02, 'T03', T03, 'T0e', T0e);

if f == 1
    vizualizare_robot(q, par);
end
end

%% Model geometric invers
function q = geometric_invers(p, par, cot)
% cot = +1 cot sus (q3 > 0, implicit), -1 cot jos (q3 < 0)
%   q1 = atan2(y, x)
%   rho = sqrt(x^2+y^2) - l3 = l4*c2 + le*c23
%   h   = l1 + l2 - z        = l4*s2 + le*s23   (mecanism plan 2R)
if nargin < 3, cot = 1; end
x = p(1); y = p(2); z = p(3);
q1  = atan2(y, x);
rho = sqrt(x^2 + y^2) - par.l3;
h   = par.l1 + par.l2 - z;
c3 = (rho^2 + h^2 - par.l4^2 - par.le^2) / (2*par.l4*par.le);
if abs(c3) > 1
    error('geometric_invers: punctul [%.3f %.3f %.3f] m este în afara spațiului de lucru.', x, y, z);
end
s3 = cot*sqrt(1 - c3^2);
q3 = atan2(s3, c3);
q2 = atan2(h, rho) - atan2(par.le*s3, par.l4 + par.le*c3);
q = [q1; q2; q3];
end

%% Jacobian
function [J, dJ] = jacobian_PUMA(q, dq, par)
% J = [Jv; Jw] (6x3):  [v; w] = J*dq,   [a; e] = J*ddq + dJ*dq
%   r = l3 + l4*c2 + le*c23,   a = l4*s2 + le*s23
%   Jv = [ -s1*r  -c1*a  -c1*le*s23 ;  c1*r  -s1*a  -s1*le*s23 ;  0  -(l4*c2+le*c23)  -le*c23 ]
%   Jw = [ 0 -s1 -s1 ;  0 c1 c1 ;  1 0 0 ]
if nargin < 2 || isempty(dq), dq = [0; 0; 0]; end
l3 = par.l3; l4 = par.l4; le = par.le;
s1 = sin(q(1)); c1 = cos(q(1));
s2 = sin(q(2)); c2 = cos(q(2));
s23 = sin(q(2)+q(3)); c23 = cos(q(2)+q(3));
dq1 = dq(1); dq2 = dq(2); dq3 = dq(3);

r  = l3 + l4*c2 + le*c23;
a  = l4*s2 + le*s23;
b  = le*s23;
cc = l4*c2 + le*c23;
c  = le*c23;

Jv = [-s1*r, -c1*a, -c1*b;
       c1*r, -s1*a, -s1*b;
       0,    -cc,   -c  ];
Jw = [0, -s1, -s1;
      0,  c1,  c1;
      1,  0,   0 ];
J = [Jv; Jw];

dr  = -a*dq2 - b*dq3;
da  =  cc*dq2 + c*dq3;
db  =  c*(dq2 + dq3);
dcc = -a*dq2 - b*dq3;
dc  = -b*(dq2 + dq3);
dJv = [-c1*dq1*r - s1*dr,   s1*dq1*a - c1*da,   s1*dq1*b - c1*db;
       -s1*dq1*r + c1*dr,  -c1*dq1*a - s1*da,  -c1*dq1*b - s1*db;
        0,                 -dcc,               -dc            ];
dJw = [0, -c1*dq1, -c1*dq1;
       0, -s1*dq1, -s1*dq1;
       0,  0,       0     ];
dJ = [dJv; dJw];
end

%% Generare traiectorie
function [pvc, t] = traiect(q, T, N, grad, flag, nume)
% q - n x 2 [start, final]; grad 3 (cubic) sau 5; pvc = [P; V; A]
%   grad 3: s = 3*tau^2 - 2*tau^3
%   grad 5: s = 10*tau^3 - 15*tau^4 + 6*tau^5     (tau = t/T)
if nargin < 4 || isempty(grad), grad = 5; end
if nargin < 5, flag = 0; end
n = size(q, 1);
if nargin < 6
    nume = arrayfun(@(k) sprintf('q_%d', k), 1:n, 'UniformOutput', false);
end
t   = linspace(0, T, N);
tau = t/T;
switch grad
    case 3
        s = 3*tau.^2 - 2*tau.^3;  ds = (6*tau - 6*tau.^2)/T;  dds = (6 - 12*tau)/T^2;
    case 5
        s   = 10*tau.^3 - 15*tau.^4 + 6*tau.^5;
        ds  = (30*tau.^2 - 60*tau.^3 + 30*tau.^4)/T;
        dds = (60*tau - 180*tau.^2 + 120*tau.^3)/T^2;
    otherwise
        error('traiect: grad trebuie să fie 3 sau 5');
end
dq = q(:,2) - q(:,1);
P = q(:,1) + dq*s;  V = dq*ds;  A = dq*dds;
pvc = [P; V; A];

if flag == 1
    figure('Name', 'Traiectorie');
    titluri = {'poziție', 'viteză', 'accelerație'};
    date = {P, V, A};
    for i = 1:n
        for k = 1:3
            subplot(n, 3, (i-1)*3 + k);
            plot(t, date{k}(i,:), 'LineWidth', 1.5); grid on;
            title([nume{i} ' - ' titluri{k}]);
            if i == n, xlabel('t [s]'); end
        end
    end
end
end

%% Model cinematic direct
function [P, V, A, W] = cinematica_directa(Q, dQ, ddQ, t, par, etapa)
% p = MGD(q),  v = Jv*dq,  a = Jv*ddq + dJv*dq,  w = Jw*dq
% etapa - numele etapei pentru grafice ('' sau lipsă = fără grafice)
if nargin < 6, etapa = ''; end
N = size(Q, 2);
P = zeros(3, N); V = P; A = P; W = P;
for i = 1:N
    P(:,i) = geometric(Q(:,i), par, 0);
    [J, dJ] = jacobian_PUMA(Q(:,i), dQ(:,i), par);
    v6 = J*dQ(:,i);
    a6 = J*ddQ(:,i) + dJ*dQ(:,i);
    V(:,i) = v6(1:3);  W(:,i) = v6(4:6);  A(:,i) = a6(1:3);
end
if ~isempty(etapa)
    tl = tab_grafice(etapa, 'Efector (cinematică directă)', 2, 2);
    g = {P, 'Poziția efectorului', {'x','y','z'}, '[m]';
         V, 'Viteza liniară', {'V_x','V_y','V_z'}, '[m/s]';
         A, 'Accelerația liniară', {'a_x','a_y','a_z'}, '[m/s^2]';
         W, 'Viteza unghiulară', {'\omega_x','\omega_y','\omega_z'}, '[rad/s]'};
    for k = 1:4
        ax = nexttile(tl);
        plot(ax, t, g{k,1}, 'LineWidth', 1.5); grid(ax, 'on');
        title(ax, g{k,2}); ylabel(ax, g{k,4});
        if k > 2, xlabel(ax, 't [s]'); end
    end
    legenda_tl(ax, {'x', 'y', 'z'}, 'componenta');
end
end

%% Model cinematic invers (urmărire r(t))
function [Q, dQ, ddQ] = cinematica_inversa(R, dR, ddR, par, cot)
% q = MGI(r),  dq = Jv^-1*dr,  ddq = Jv^-1*(ddr - dJv*dq)
if nargin < 5, cot = 1; end
N = size(R, 2);
Q = zeros(3, N); dQ = Q; ddQ = Q;
for i = 1:N
    q = geometric_invers(R(:,i), par, cot);
    J = jacobian_PUMA(q, [], par);
    Jv = J(1:3,:);
    if rcond(Jv) < 1e-8
        warning('cinematica_inversa: configurație singulară la pasul %d', i);
    end
    dq = Jv \ dR(:,i);
    [~, dJ] = jacobian_PUMA(q, dq, par);
    ddq = Jv \ (ddR(:,i) - dJ(1:3,:)*dq);
    Q(:,i) = q; dQ(:,i) = dq; ddQ(:,i) = ddq;
end
Q(1,:) = unwrap(Q(1,:));   % evită saltul de 2*pi al lui atan2
end

%% Model dinamic
function [M, C, G, V] = model_dinamic(q, dq, par)
% M(q)*ddq + C(q,dq)*dq + G(q) + Fv*dq = tau     (Lagrange)
%   r_c2 = l3 + lc2*c2,           z_c2 = l1 + l2 - lc2*s2
%   r_c3 = l3 + l4*c2 + lc3*c23,  z_c3 = l1 + l2 - l4*s2 - lc3*s23
%   w2 = [-s2*dq1; dq2; c2*dq1],  w3 = [-s23*dq1; dq2+dq3; c23*dq1]  (în axele proprii)
if nargin < 2 || isempty(dq), dq = [0; 0; 0]; end
l1 = par.l1; l2 = par.l2; l3 = par.l3; l4 = par.l4;
m2 = par.m2; m3 = par.m3; lc2 = par.lc2; lc3 = par.lc3; g = par.g;
I2xx = par.I2(1); I2yy = par.I2(2); I2zz = par.I2(3);
I3xx = par.I3(1); I3yy = par.I3(2); I3zz = par.I3(3);
s2 = sin(q(2)); c2 = cos(q(2)); s3 = sin(q(3)); c3 = cos(q(3));
s23 = sin(q(2)+q(3)); c23 = cos(q(2)+q(3));
dq1 = dq(1); dq2 = dq(2); dq3 = dq(3);

rc2 = l3 + lc2*c2;
rc3 = l3 + l4*c2 + lc3*c23;

% Matricea maselor
M11 = par.I1zz + I2xx*s2^2 + I2zz*c2^2 + I3xx*s23^2 + I3zz*c23^2 + m2*rc2^2 + m3*rc3^2;
M22 = I2yy + m2*lc2^2 + I3yy + m3*(l4^2 + lc3^2 + 2*l4*lc3*c3);
M23 = I3yy + m3*(lc3^2 + l4*lc3*c3);
M33 = I3yy + m3*lc3^2;
M = [M11 0   0;
     0   M22 M23;
     0   M23 M33];

% Matricea Coriolis/centrifugală (Christoffel)
% h2 = 1/2*dM11/dq2,  h3 = 1/2*dM11/dq3,  h = m3*l4*lc3*s3
h2 = (I2xx - I2zz)*s2*c2 + (I3xx - I3zz)*s23*c23 - m2*rc2*lc2*s2 - m3*rc3*(l4*s2 + lc3*s23);
h3 = (I3xx - I3zz)*s23*c23 - m3*rc3*lc3*s23;
h  = m3*l4*lc3*s3;
C = [ h2*dq2 + h3*dq3,  h2*dq1,   h3*dq1;
     -h2*dq1,          -h*dq3,   -h*(dq2 + dq3);
     -h3*dq1,           h*dq2,    0            ];

% Vectorul gravitațional G = dV/dq
G = [0;
     -g*(m2*lc2*c2 + m3*(l4*c2 + lc3*c23));
     -g*m3*lc3*c23];

% Energia potențială
V = g*(m2*(l1 + l2 - lc2*s2) + m3*(l1 + l2 - l4*s2 - lc3*s23));
end

%% Dinamica inversă
function Tau = dinamica_inversa(Q, dQ, ddQ, t, par, etapa)
% tau = M*ddq + C*dq + G + Fv*dq
% etapa - numele etapei pentru grafice ('' sau lipsă = fără grafice)
if nargin < 6, etapa = ''; end
N = size(Q, 2);
Tau = zeros(3, N); Tin = Tau; Tcor = Tau; Tgr = Tau;
for i = 1:N
    [M, C, G] = model_dinamic(Q(:,i), dQ(:,i), par);
    Tin(:,i)  = M*ddQ(:,i);
    Tcor(:,i) = C*dQ(:,i);
    Tgr(:,i)  = G;
    Tau(:,i)  = Tin(:,i) + Tcor(:,i) + Tgr(:,i) + par.Fv(:).*dQ(:,i);
end
if ~isempty(etapa)
    tl = tab_grafice(etapa, 'Cupluri (dinamică inversă)', 3, 1);
    for k = 1:3
        ax = nexttile(tl); hold(ax, 'on');
        plot(ax, t, Tau(k,:), 'k', 'LineWidth', 2);
        plot(ax, t, Tin(k,:), '--', t, Tcor(k,:), '--', t, Tgr(k,:), '--', 'LineWidth', 1.2);
        grid(ax, 'on'); ylabel(ax, sprintf('\\tau_%d [N·m]', k));
        if k == 1, title(ax, 'Cuplurile din cuple'); end
    end
    xlabel(ax, 't [s]');
    legenda_tl(ax, {'\tau total', 'M·ddq', 'C·dq', 'G'});
end
end

%% Dinamica directă
function dx = dinamica_directa(t, x, tau_fun, par)
% x = [q; dq],  ddq = M^-1*(tau - C*dq - G - Fv*dq)
% tau_fun(t, x): cuplu în buclă deschisă, zero, sau (mai târziu) legea de reglare
q  = x(1:3);
dq = x(4:6);
tau = tau_fun(t, x);
[M, C, G] = model_dinamic(q, dq, par);
ddq = M \ (tau(:) - C*dq - G - par.Fv(:).*dq);
dx = [dq; ddq];
end

function [Q, dQ] = simulare_dinamica(x0, t, tau_fun, par)
% integrarea modelului dinamic cu ode45
opt = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[~, X] = ode45(@(tt, x) dinamica_directa(tt, x, tau_fun, par), t, x0(:), opt);
Q  = X(:, 1:3).';
dQ = X(:, 4:6).';
end

%% ===================== CONTROL PI (ETAPA 5) ============================
function [t, Q, R] = referinta_control(date, par)
% Traiectoria dorită a cuplelor q_d(t), eșantionată cu pasul dt_control:
%   'r' - din traiectoria efectorului r(t) prin r_puncte (MGI + cinematică inversă)
%   'q' - direct prin punctele q_puncte (spațiul cuplelor)
R = [];
if strcmpi(date.control_ref, 'r') && size(date.r_puncte, 1) >= 2
    [t, rvc] = traiectorie_puncte(date.r_puncte.', date.T_segment, [], date.grad, date.dt_control);
    R = rvc(1:3,:);
    verifica_spatiu_lucru(R, t, par);
    Q = cinematica_inversa(R, rvc(4:6,:), rvc(7:9,:), par, date.cot);
    fprintf('Referința: traiectoria r(t) prin %d puncte, %.2f s (q_d din modelul cinematic invers)\n', ...
            size(date.r_puncte, 1), t(end));
elseif size(date.q_puncte, 1) >= 2
    [t, pvc] = traiectorie_puncte(deg2rad(date.q_puncte).', date.T_segment, [], date.grad, date.dt_control);
    Q = pvc(1:3,:);
    fprintf('Referința: mișcare prin %d puncte în spațiul cuplelor, %.2f s\n', size(date.q_puncte, 1), t(end));
else
    error('Pentru control trebuie date r_puncte (control_ref = ''r'') sau q_puncte (control_ref = ''q'').');
end
end

function act = actionari(date, par, Q)
% Parametrii acționărilor (motor cc + reductor) și acordarea regulatoarelor PI.
% Modelul motorului (ca la Pozna), pentru fiecare cuplă:
%   ddq = k2*v - k1*dq - k3*tau
%   k1 = (bm + km^2/R)/Ja,   k2 = km/(Ja*ir*R),   k3 = 1/(Ja*ran*ir^2)
% Funcția de transfer folosită la acordare (tau tratat ca perturbație):
%   G(s) = k2/(s*(s + k1)),   C(s) = (Kp*s + Ki)/s
%   P(s) = s^3 + k1*s^2 + k2*Kp*s + k2*Ki = (s + p1)(s + p2)(s + p3)
%   => p3 = k1 - p1 - p2,  Kp = (p1*p2 + p1*p3 + p2*p3)/k2,  Ki = p1*p2*p3/k2
% În k1 și k2 de la acordare se folosește inerția echivalentă (motor + braț):
%   J_ech = Ja + M_ii/(ran*ir^2)
% cu M_ii (elementul diagonal al matricei maselor) mediu ('medie') sau maxim
% ('max') pe traiectorie; cu 'motor' se neglijează brațul (J_ech = Ja).
m = date.motor;
c = @(v) v(:);                                     % coloană (câte o valoare pe cuplă)
act.nume = m.nume;
act.U_nom = c(m.U_nom);  act.R = c(m.R);  act.km = c(m.km);  act.ir = c(m.ir);
act.ran = c(m.ran);  act.M_nom = c(m.M_nom);  act.I_nom = c(m.I_nom);
act.M_red = c(m.M_red);  act.n_red = c(m.n_red);  act.n0 = c(m.n0);
act.Ja = c(m.J_rotor) + c(m.J_red);                % rotor + reductor (la arborele motorului)
act.bm = act.km.*c(m.I0)./(c(m.n0)*pi/30);         % frecare vâscoasă estimată din mersul în gol
act.k1 = (act.bm + act.km.^2./act.R)./act.Ja;
act.k2 = act.km./(act.Ja.*act.ir.*act.R);
act.k3 = 1./(act.Ja.*act.ran.*act.ir.^2);

% Inerția brațului văzută de fiecare cuplă: M_ii(q) pe toată traiectoria
Mt = zeros(3, size(Q, 2));
for i = 1:size(Q, 2)
    Mt(:,i) = diag(model_dinamic(Q(:,i), [], par));
end
act.Mii = [min(Mt, [], 2), mean(Mt, 2), max(Mt, [], 2)];   % [min medie max]
act.raport_inertii = act.Mii.*act.k3;                       % k3*M_ii (la Pozna ~0.05)
act.acordare = lower(date.acordare);
switch act.acordare
    case 'medie', Md = act.Mii(:,2);
    case 'max',   Md = act.Mii(:,3);
    case 'motor', Md = zeros(3, 1);
    otherwise, error('acordare trebuie să fie ''medie'', ''max'' sau ''motor''.');
end
act.J_ech = act.Ja + Md./(act.ran.*act.ir.^2);
act.k1e = (act.bm + act.km.^2./act.R)./act.J_ech;
act.k2e = act.km./(act.J_ech.*act.ir.*act.R);

% Alocarea polilor (poli reali în -p1, -p2, -p3)
p1 = date.poli(1)*act.k1e;  p2 = date.poli(2)*act.k1e;  p3 = act.k1e - p1 - p2;
act.poli = [p1 p2 p3];
act.Kp = (p1.*p2 + p1.*p3 + p2.*p3)./act.k2e;
act.Ki = p1.*p2.*p3./act.k2e;

% Tensiunea care menține brațul în poziția inițială (starea inițială a integratorului PI)
[~, ~, G0] = model_dinamic(Q(:,1), [], par);
act.v0 = act.R.*G0./(act.km.*act.ir.*act.ran);

% Polii reali ai buclei (cu inerția brațului M_ii(q) de pe traiectorie):
%   (1 + k3*M_ii)*ddq = k2*v - k1*dq - ...  =>  k1 și k2 se împart la 1 + k3*M_ii
act.zeta_min = ones(3, 1);
for i = 1:3
    for Mi = linspace(act.Mii(i,1), act.Mii(i,3), 11)
        a = 1 + act.k3(i)*Mi;
        r = roots([1, act.k1(i)/a, act.k2(i)*act.Kp(i)/a, act.k2(i)*act.Ki(i)/a]);
        act.zeta_min(i) = min([act.zeta_min(i); -real(r)./abs(r)]);
    end
end
end

function afiseaza_actionari(act)
fprintf('Acționări (motor cc + reductor) și regulatoare PI (alocarea polilor, acordare ''%s''):\n', act.acordare);
for i = 1:3
    fprintf(' Cupla %d: %s\n', i, act.nume{i});
    fprintf('   Ja = %.3e kg·m^2, bm = %.2e N·m·s/rad, ir = %g, ran = %.2f\n', ...
            act.Ja(i), act.bm(i), act.ir(i), act.ran(i));
    fprintf('   k1 = %.2f 1/s, k2 = %.3f rad/(V·s^2), k3 = %.3f 1/(kg·m^2)\n', act.k1(i), act.k2(i), act.k3(i));
    fprintf('   M_ii = %.3f ... %.3f kg·m^2 (medie %.3f), k3·M_ii = %.2f ... %.2f, J_ech = %.3e kg·m^2\n', ...
            act.Mii(i,1), act.Mii(i,3), act.Mii(i,2), act.raport_inertii(i,1), act.raport_inertii(i,3), act.J_ech(i));
    fprintf('   acordare: k1 = %.2f, k2 = %.3f -> poli = -[%.2f %.2f %.2f] rad/s\n', ...
            act.k1e(i), act.k2e(i), act.poli(i,:));
    fprintf('   Kp = %.4g V/rad,  Ki = %.4g V/(rad·s),  v0 = %.3f V,  amortizare reală min = %.2f\n', ...
            act.Kp(i), act.Ki(i), act.v0(i), act.zeta_min(i));
    r = roots([1, act.k1e(i), act.k2e(i)*act.Kp(i), act.k2e(i)*act.Ki(i)]);
    if max(abs(sort(real(r)) + sort(act.poli(i,:), 'descend').')) > 1e-3*act.k1e(i)
        warning('Polinomul caracteristic nu are rădăcinile alese (cupla %d).', i);
    end
end
end

function res = control_matlab(tc, Qc, act, par, t_stop)
% Bucla de reglare simulată direct în MATLAB (verificare pentru modelul Simulink).
% Cuplajul motor-structură se rezolvă exact (fără bucla algebrică):
%   ddq = k2.*v - k1.*dq - k3.*tau,   tau = M*ddq + C*dq + G + Fv.*dq
%   => (I + diag(k3)*M)*ddq = k2.*v - k1.*dq - k3.*(C*dq + G + Fv.*dq)
Fq = griddedInterpolant(tc, Qc.', 'linear', 'nearest');
x0 = [Qc(:,1); 0; 0; 0; act.v0];          % [q; dq; xi], xi = starea integratorului PI [V]
t  = 0:0.005:t_stop;
opt = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 1e-2);
[~, X] = ode45(@(tt, x) bucla_control(tt, x, Fq, act, par), t, x0, opt);
X = X.';
res.t = t;  res.Q = X(1:3,:);  res.dQ = X(4:6,:);
res.ddQ = zeros(3, numel(t));  res.Tau = res.ddQ;  res.U = res.ddQ;  res.E = res.ddQ;
for i = 1:numel(t)
    [~, res.ddQ(:,i), res.Tau(:,i), res.U(:,i), res.E(:,i)] = bucla_control(t(i), X(:,i), Fq, act, par);
end
end

function [dx, ddq, tau, v, e] = bucla_control(t, x, Fq, act, par)
q = x(1:3);  dq = x(4:6);  xi = x(7:9);
e = Fq(t).' - q;                                            % eroarea de poziție
v = min(max(act.Kp.*e + xi, -act.U_nom), act.U_nom);        % PI + saturație
[M, C, G] = model_dinamic(q, dq, par);
h = C*dq + G + par.Fv.*dq;
ddq = (eye(3) + diag(act.k3)*M) \ (act.k2.*v - act.k1.*dq - act.k3.*h);
tau = M*ddq + h;
dx = [dq; ddq; act.Ki.*e];
end

function res = control_simulink(fisier, act, par, q0, t_stop, date, folder)
% Construiește modelul Simulink, îl salvează și îl simulează.
% Modelul se regenerează la fiecare rulare (din datele din Date_PUMA.m).
[~, mdl] = fileparts(fisier);
dir0 = cd(folder);                 % From File citește Poz_dorita.mat din folderul proiectului
revino = onCleanup(@() cd(dir0));
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist(fisier, 'file'), delete(fisier); end
construieste_model(mdl, act, par, q0, t_stop, date.rupere_bucla, date.T_filtru);
save_system(mdl, fisier);
out = sim(mdl);
if date.deschide_simulink && ~batchStartupOptionUsed
    open_system(mdl);              % modelul se deschide în Simulink, ca la Pozna
end
[res.t, iu] = unique(out.miscari.Time(:).', 'last');     % fără momente de timp duplicate
serie = @(s) subsref(reshape(s.Data, numel(s.Time), []).', struct('type', '()', 'subs', {{':', iu}}));
mis = serie(out.miscari);
res.Q   = mis([1 4 7],:);  res.dQ = mis([2 5 8],:);  res.ddQ = mis([3 6 9],:);
res.Tau = serie(out.moMo);
res.E   = serie(out.erori);
res.U   = serie(out.comanda);
end

function construieste_model(mdl, act, par, q0, t_stop, rupere, Tf)
% Modelul Simulink, cu aceeași structură ca modelul contrClasic42 (prof. Pozna):
%   Poz_dorita (From File) -> Controler -> MotoReductor1..3 -> Bus -> StructuraPUMA -> Tau
% Bucla algebrică ddq -> Tau -> ddq (motorul are nevoie de Tau, structura de ddq)
% se rupe în StructuraPUMA:
%   'filtru' - filtru de ordinul 1, Tau_f' = (Tau - Tau_f)/Tf, cu Tf mic (stabil oricare
%              ar fi raportul de inerții; sugestia Simulink pentru bucle continue)
%   'memory' - blocul Memory, ca la Pozna (stabil doar dacă k3*M_ii < 1)
new_system(mdl);
set_param(mdl, 'StopTime', num2str(t_stop), 'Solver', 'ode45', 'RelTol', '1e-6', ...
          'AbsTol', '1e-8', 'MaxStep', '1e-2', 'AlgebraicLoopMsg', 'none', ...
          'ReturnWorkspaceOutputs', 'on');
num = @(v) sprintf('%.10g', v);
vec = @(v) ['[' strjoin(arrayfun(num, v(:).', 'UniformOutput', false), ' ') ']'];
blk = @(lib, cale, poz, varargin) add_block(lib, cale, 'Position', poz, varargin{:});
lin = @(sis, a, b) add_line(sis, a, b, 'autorouting', 'smart');

% ---- nivelul principal
blk('simulink/Sources/From File', [mdl '/Poz_dorita'], [30 300 150 360], ...
    'FileName', 'Poz_dorita.mat', 'ExtrapolationBeforeFirstDataPoint', 'Hold first value', ...
    'ExtrapolationAfterLastDataPoint', 'Hold last value');
blk('built-in/Subsystem', [mdl '/Controler'], [240 120 370 520]);
for i = 1:3
    blk('built-in/Subsystem', sprintf('%s/MotoReductor%d', mdl, i), [490 60+180*(i-1) 620 160+180*(i-1)]);
end
blk('simulink/Signal Routing/Mux', [mdl '/Bus'], [720 40 728 520], 'Inputs', '9', 'DisplayOption', 'bar');
blk('simulink/Sinks/To Workspace', [mdl '/out.miscari'], [800 30 890 60], ...
    'VariableName', 'miscari', 'SaveFormat', 'Timeseries');
blk('built-in/Subsystem', [mdl '/StructuraPUMA'], [800 240 940 320]);
blk('simulink/Signal Routing/Demux', [mdl '/Tau'], [1000 240 1008 320], 'Outputs', '3', 'DisplayOption', 'bar');
blk('simulink/Sinks/To Workspace', [mdl '/out.moMo'], [1050 380 1140 410], ...
    'VariableName', 'moMo', 'SaveFormat', 'Timeseries');
blk('simulink/Signal Routing/Mux', [mdl '/Mux_n'], [720 580 728 660], 'Inputs', '3', 'DisplayOption', 'bar');
blk('simulink/Sinks/To Workspace', [mdl '/out.turatii'], [800 605 890 635], ...
    'VariableName', 'turatii', 'SaveFormat', 'Timeseries');

% ---- Controler: e = q_d - q -> PI -> saturație -> v
s = [mdl '/Controler'];
for i = 1:3
    blk('built-in/Inport', sprintf('%s/q%d', s, i), [30 40+150*(i-1) 60 54+150*(i-1)], 'Port', num2str(i), ...
        'PortDimensions', '1');
end
blk('built-in/Inport', [s '/q1d_q2d_q3d'], [30 480 60 494], 'Port', '4', 'PortDimensions', '3');
blk('simulink/Signal Routing/Demux', [s '/Demux'], [120 440 125 540], 'Outputs', '3', 'DisplayOption', 'bar');
blk('simulink/Signal Routing/Mux', [s '/Mux_e'], [400 520 405 600], 'Inputs', '3', 'DisplayOption', 'bar');
blk('simulink/Signal Routing/Mux', [s '/Mux_v'], [640 520 645 600], 'Inputs', '3', 'DisplayOption', 'bar');
blk('simulink/Sinks/To Workspace', [s '/out.erori'], [450 545 540 575], 'VariableName', 'erori', 'SaveFormat', 'Timeseries');
blk('simulink/Sinks/To Workspace', [s '/out.comanda'], [690 545 780 575], 'VariableName', 'comanda', 'SaveFormat', 'Timeseries');
lin(s, 'q1d_q2d_q3d/1', 'Demux/1');
for i = 1:3
    y = 40 + 150*(i-1);
    blk('simulink/Math Operations/Sum', sprintf('%s/e%d', s, i), [200 y-5 225 y+20], 'Inputs', '+-', 'IconShape', 'round');
    blk('simulink/Continuous/PID Controller', sprintf('%s/PI%d', s, i), [280 y-10 350 y+25], ...
        'Controller', 'PI', 'P', num(act.Kp(i)), 'I', num(act.Ki(i)), ...
        'InitialConditionForIntegrator', num(act.v0(i)));
    blk('simulink/Discontinuities/Saturation', sprintf('%s/Saturatie%d', s, i), [400 y-10 440 y+25], ...
        'UpperLimit', num(act.U_nom(i)), 'LowerLimit', num(-act.U_nom(i)));
    blk('built-in/Outport', sprintf('%s/v%d', s, i), [560 y 590 y+14], 'Port', num2str(i));
    lin(s, sprintf('Demux/%d', i), sprintf('e%d/1', i));
    lin(s, sprintf('q%d/1', i), sprintf('e%d/2', i));
    lin(s, sprintf('e%d/1', i), sprintf('PI%d/1', i));
    lin(s, sprintf('e%d/1', i), sprintf('Mux_e/%d', i));
    lin(s, sprintf('PI%d/1', i), sprintf('Saturatie%d/1', i));
    lin(s, sprintf('Saturatie%d/1', i), sprintf('v%d/1', i));
    lin(s, sprintf('Saturatie%d/1', i), sprintf('Mux_v/%d', i));
end
lin(s, 'Mux_e/1', 'out.erori/1');
lin(s, 'Mux_v/1', 'out.comanda/1');

% ---- MotoReductor i: ddq = k2*v - k1*dq - k3*tau -> 1/s -> dq -> 1/s -> q
for i = 1:3
    s = sprintf('%s/MotoReductor%d', mdl, i);
    set_param(s, 'AttributesFormatString', act.nume{i});
    blk('built-in/Inport', [s '/Tau'], [30 100 60 114], 'Port', '1', 'PortDimensions', '1');
    blk('built-in/Inport', [s '/u'], [30 220 60 234], 'Port', '2', 'PortDimensions', '1');
    blk('simulink/Signal Routing/Mux', [s '/Mux'], [130 80 135 240], 'Inputs', '3', 'DisplayOption', 'bar');
    blk('simulink/User-Defined Functions/MATLAB Function', [s '/Acc MotoRed'], [190 130 300 190]);
    seteaza_cod([s '/Acc MotoRed'], cod_motor(act, i), {'U', '3'; 'ddq', '1'});
    blk('simulink/Continuous/Integrator', [s '/Integrator_dq'], [360 145 390 175], 'InitialCondition', '0');
    blk('simulink/Continuous/Integrator', [s '/Integrator_q'], [460 145 490 175], 'InitialCondition', num(q0(i)));
    blk('simulink/Math Operations/Gain', [s '/rps2rpm'], [460 230 510 260], 'Gain', num(act.ir(i)*30/pi));
    blk('built-in/Outport', [s '/q'], [570 153 600 167], 'Port', '1');
    blk('built-in/Outport', [s '/dq'], [570 60 600 74], 'Port', '2');
    blk('built-in/Outport', [s '/ddq'], [570 320 600 334], 'Port', '3');
    blk('built-in/Outport', [s '/n'], [570 238 600 252], 'Port', '4');
    lin(s, 'Integrator_dq/1', 'Mux/1');
    lin(s, 'Tau/1', 'Mux/2');
    lin(s, 'u/1', 'Mux/3');
    lin(s, 'Mux/1', 'Acc MotoRed/1');
    lin(s, 'Acc MotoRed/1', 'Integrator_dq/1');
    lin(s, 'Acc MotoRed/1', 'ddq/1');
    lin(s, 'Integrator_dq/1', 'Integrator_q/1');
    lin(s, 'Integrator_dq/1', 'dq/1');
    lin(s, 'Integrator_dq/1', 'rps2rpm/1');
    lin(s, 'Integrator_q/1', 'q/1');
    lin(s, 'rps2rpm/1', 'n/1');
end

% ---- StructuraPUMA: modelul dinamic invers (q, dq, ddq -> Tau)
s = [mdl '/StructuraPUMA'];
blk('built-in/Inport', [s '/u'], [30 98 60 112], 'Port', '1', 'PortDimensions', '9');
blk('simulink/User-Defined Functions/MATLAB Function', [s '/Dinamica inversa PUMA'], [200 75 340 135]);
seteaza_cod([s '/Dinamica inversa PUMA'], cod_structura(par), {'u', '9'; 'Tau', '3'});
blk('built-in/Outport', [s '/Tau'], [700 98 730 112], 'Port', '1');
switch lower(rupere)
    case 'memory'
        % ca la Pozna: blocul Memory (întârziere de un pas de integrare)
        u0 = [q0(1) 0 0 q0(2) 0 0 q0(3) 0 0];
        blk('simulink/Discrete/Memory', [s '/Memory'], [110 90 140 120], 'InitialCondition', vec(u0));
        lin(s, 'u/1', 'Memory/1');
        lin(s, 'Memory/1', 'Dinamica inversa PUMA/1');
        lin(s, 'Dinamica inversa PUMA/1', 'Tau/1');
    case 'filtru'
        % Tau_f = Tau/(Tf*s + 1), pornește din cuplul static din poziția inițială
        [~, ~, G0] = model_dinamic(q0, [], par);
        lin(s, 'u/1', 'Dinamica inversa PUMA/1');
        blk('simulink/Math Operations/Sum', [s '/Sum_f'], [400 95 420 115], 'Inputs', '+-', 'IconShape', 'round');
        blk('simulink/Math Operations/Gain', [s '/1_Tf'], [460 90 510 120], 'Gain', num(1/Tf));
        blk('simulink/Continuous/Integrator', [s '/Filtru Tau'], [560 90 590 120], 'InitialCondition', vec(G0));
        set_param([s '/Filtru Tau'], 'AttributesFormatString', sprintf('Tf = %g s', Tf));
        lin(s, 'Dinamica inversa PUMA/1', 'Sum_f/1');
        lin(s, 'Sum_f/1', '1_Tf/1');
        lin(s, '1_Tf/1', 'Filtru Tau/1');
        lin(s, 'Filtru Tau/1', 'Sum_f/2');
        lin(s, 'Filtru Tau/1', 'Tau/1');
    otherwise
        error('rupere_bucla trebuie să fie ''filtru'' sau ''memory''.');
end

% ---- legăturile de la nivelul principal
lin(mdl, 'Poz_dorita/1', 'Controler/4');
for i = 1:3
    m = sprintf('MotoReductor%d', i);
    lin(mdl, sprintf('Controler/%d', i), [m '/2']);
    lin(mdl, [m '/1'], sprintf('Controler/%d', i));
    lin(mdl, [m '/1'], sprintf('Bus/%d', 3*i - 2));
    lin(mdl, [m '/2'], sprintf('Bus/%d', 3*i - 1));
    lin(mdl, [m '/3'], sprintf('Bus/%d', 3*i));
    lin(mdl, [m '/4'], sprintf('Mux_n/%d', i));
    lin(mdl, sprintf('Tau/%d', i), [m '/1']);
end
lin(mdl, 'Bus/1', 'StructuraPUMA/1');
lin(mdl, 'Bus/1', 'out.miscari/1');
lin(mdl, 'StructuraPUMA/1', 'Tau/1');
lin(mdl, 'StructuraPUMA/1', 'out.moMo/1');
lin(mdl, 'Mux_n/1', 'out.turatii/1');
end

function seteaza_cod(bloc, cod, dimensiuni)
% Scrie codul în blocul MATLAB Function și fixează dimensiunile intrărilor/ieșirilor
% (dimensiuni = {nume, 'n'; ...}), ca Simulink să nu le mai deducă prin bucla de reacție.
rt = sfroot;
ch = rt.find('-isa', 'Stateflow.EMChart', 'Path', bloc);
ch.Script = cod;
for k = 1:size(dimensiuni, 1)
    d = ch.find('-isa', 'Stateflow.Data', 'Name', dimensiuni{k, 1});
    if isempty(d), error('Blocul %s nu are variabila %s.', bloc, dimensiuni{k, 1}); end
    d.Props.Array.Size = dimensiuni{k, 2};
end
end

function cod = cod_motor(act, i)
% Codul blocului "Acc MotoRed" (același model ca la Pozna)
cod = strjoin({
    'function ddq = acc_motored(U)'
    sprintf('%% Motor cc + reductor, cupla %d: %s', i, act.nume{i})
    '% U = [dq; tau; v]: viteza cuplei, cuplul cerut de structură, tensiunea'
    'dq  = U(1);'
    'tau = U(2);'
    'v   = U(3);'
    '% Parametrii motor + reductor (din Date_PUMA.m)'
    sprintf('Ja  = %.10g;   %% inerția rotor + reductor [kg*m^2]', act.Ja(i))
    sprintf('R   = %.10g;   %% rezistența bobinajului [Ohm]', act.R(i))
    sprintf('bm  = %.10g;   %% frecarea vâscoasă din motor [N*m*s/rad]', act.bm(i))
    sprintf('km  = %.10g;   %% constanta de cuplu [N*m/A]', act.km(i))
    sprintf('ir  = %.10g;   %% raportul de transmisie', act.ir(i))
    sprintf('ran = %.10g;   %% randamentul reductorului', act.ran(i))
    '% Coeficienții'
    'k1 = (bm + km^2/R)/Ja;   % dq'
    'k2 = km/(Ja*ir*R);       % v'
    'k3 = 1/(Ja*ran*ir^2);    % tau'
    '% Ecuația'
    'ddq = k2*v - k1*dq - k3*tau;'
    'end'}, newline);
end

function cod = cod_structura(par)
% Codul blocului "Dinamica inversa PUMA": folosește exact funcția model_dinamic
% din acest fișier (copiată în bloc), deci modelul din Simulink este același.
sursa = fileread([mfilename('fullpath') '.m']);
md = regexp(sursa, 'function \[M, C, G, V\] = model_dinamic.*?\nend(?=\s*\n)', 'match', 'once');
if isempty(md), error('Nu am găsit funcția model_dinamic în %s.m', mfilename); end
num = @(v) sprintf('%.10g', v);
vec = @(v) ['[' strjoin(arrayfun(num, v(:).', 'UniformOutput', false), ' ') ']'];
cod = strjoin({
    'function Tau = dinamica_inversa_PUMA(u)'
    '% Modelul dinamic invers al structurii PUMA (Lagrange):'
    '%   Tau = M(q)*ddq + C(q,dq)*dq + G(q) + Fv.*dq'
    '% u = [q1 dq1 ddq1 q2 dq2 ddq2 q3 dq3 ddq3]'
    'q   = [u(1); u(4); u(7)];'
    'dq  = [u(2); u(5); u(8)];'
    'ddq = [u(3); u(6); u(9)];'
    '% Parametrii structurii (din Date_PUMA.m)'
    ['par.l1 = ' num(par.l1) ';  par.l2 = ' num(par.l2) ';  par.l3 = ' num(par.l3) ...
     ';  par.l4 = ' num(par.l4) ';  par.g = ' num(par.g) ';']
    ['par.m2 = ' num(par.m2) ';  par.m3 = ' num(par.m3) ';  par.lc2 = ' num(par.lc2) ...
     ';  par.lc3 = ' num(par.lc3) ';']
    ['par.I1zz = ' num(par.I1zz) ';']
    ['par.I2 = ' vec(par.I2) ';   % [Ixx Iyy Izz] element 2']
    ['par.I3 = ' vec(par.I3) ';   % [Ixx Iyy Izz] element 3 + sarcina']
    ['Fv = ' vec(par.Fv) ''';  % frecarea vâscoasă în cuple']
    '[M, C, G] = model_dinamic(q, dq, par);'
    'Tau = M*ddq + C*dq + G + Fv.*dq;'
    'end'
    ''
    md}, newline);
end

function res = marimi_motor(res, tc, Qc, act)
% Mărimile din motoare, calculate din mișcarea simulată
res.Qd = interp1(tc, Qc.', min(res.t, tc(end))).';
res.n  = act.ir.*res.dQ*30/pi;                                % turația motoarelor [rpm]
res.I  = (res.U - act.km.*act.ir.*res.dQ)./act.R;             % curentul (fără inductanță) [A]
res.Mm = act.km.*res.I;                                       % cuplul motorului [N·m]
end

function verifica_actionari(res, act)
% Verificarea dimensionării pe mișcarea simulată (valori max / efective)
rms_ = @(x) sqrt(trapz(res.t, x.^2, 2)/(res.t(end) - res.t(1)));
fprintf('Verificarea acționărilor pe mișcarea simulată:\n');
fprintf('  cupla | max|e| [rad] | max|v| / U_nom [V] | I_rms / I_nom [A] | n_max / n_red [rpm] | max|tau| / M_red [N·m]\n');
for i = 1:3
    ok = max(abs(res.U(i,:))) < act.U_nom(i) && rms_(res.I(i,:)) <= act.I_nom(i) && ...
         max(abs(res.n(i,:))) <= act.n_red(i) && max(abs(res.Tau(i,:))) <= act.M_red(i);
    stare = 'OK';  if ~ok, stare = 'DEPĂȘIT'; end
    fprintf('    %d   |   %.2e   |    %5.2f / %4.1f     |   %5.2f / %5.2f   |   %6.0f / %5.0f    |  %6.2f / %5.1f   %s\n', ...
            i, max(abs(res.E(i,:))), max(abs(res.U(i,:))), act.U_nom(i), rms_(res.I(i,:)), act.I_nom(i), ...
            max(abs(res.n(i,:))), act.n_red(i), max(abs(res.Tau(i,:))), act.M_red(i), stare);
end
end

function grafice_control(res, tc, Qc, act, sursa)
% Graficele controlului (ca la Pozna: q, dq, ddq, tau, e, u), în tabul "5. Control PI"
et = '5. Control PI';
t = res.t;
cul = lines(3);
cu3 = @(ax, Y, lw) linii3(ax, t, Y, cul, lw);
leg = @(ax, s) legenda_tl(ax, arrayfun(@(k) sprintf('%s%d', s, k), 1:3, 'UniformOutput', false));

% 1. Pozițiile: obținute (linie plină) și dorite (punctat)
ax = nexttile(tab_grafice(et, 'Poziții', 1, 1)); hold(ax, 'on');
cu3(ax, res.Q, 2);
for k = 1:3, plot(ax, tc, Qc(k,:), '--k', 'LineWidth', 1); end
grid(ax, 'on'); xlabel(ax, 'timp [s]'); ylabel(ax, 'q1, q2, q3 [rad]');
legenda_tl(ax, {'q1', 'q2', 'q3', 'q_d (dorit)'});
title(ax, ['Poziția cuplelor - ' sursa]);

% 2. Vitezele și accelerațiile
tl = tab_grafice(et, 'Viteze și accelerații', 2, 1);
ax = nexttile(tl); hold(ax, 'on'); cu3(ax, res.dQ, 1.5);
grid(ax, 'on'); ylabel(ax, 'dq1, dq2, dq3 [rad/s]'); title(ax, 'Viteza cuplelor');
ax = nexttile(tl); hold(ax, 'on'); cu3(ax, res.ddQ, 1.5);
grid(ax, 'on'); ylabel(ax, 'ddq1, ddq2, ddq3 [rad/s^2]'); xlabel(ax, 'timp [s]');
leg(ax, 'cupla ');
title(ax, 'Accelerația cuplelor');

% 3-5. Momentele, erorile, comanda
g = {res.Tau, 'Momente', 'tau', 'tau1, tau2, tau3 [Nm]', 'Momentele din cuple (StructuraPUMA)';
     res.E,   'Erori',   'e',   'e1, e2, e3 [rad]',      'Erorile de poziție (e = q_d - q)';
     res.U,   'Comanda', 'u',   'u1, u2, u3 [V]',        'Comanda (tensiunea motoarelor)'};
for j = 1:size(g, 1)
    ax = nexttile(tab_grafice(et, g{j, 2}, 1, 1)); hold(ax, 'on');
    cu3(ax, g{j, 1}, 1.5);
    grid(ax, 'on'); xlabel(ax, 'timp [s]'); ylabel(ax, g{j, 4}); leg(ax, g{j, 3}); title(ax, g{j, 5});
end

% 6. Motoarele: curent, cuplu, turație
tl = tab_grafice(et, 'Motoare', 3, 1);
ax = nexttile(tl); hold(ax, 'on'); cu3(ax, res.I, 1.5);
grid(ax, 'on'); ylabel(ax, 'I [A]'); title(ax, 'Curentul, cuplul și turația motoarelor');
ax = nexttile(tl); hold(ax, 'on'); cu3(ax, 1000*res.Mm, 1.5);
plot(ax, t([1 end]), 1000*act.M_nom(1)*[1 1], ':k', t([1 end]), -1000*act.M_nom(1)*[1 1], ':k');
grid(ax, 'on'); ylabel(ax, 'M_{motor} [mNm]');
legenda_tl(ax, {'motor 1', 'motor 2', 'motor 3', '\pmM_{nominal}'});
ax = nexttile(tl); hold(ax, 'on'); cu3(ax, res.n, 1.5);
grid(ax, 'on'); ylabel(ax, 'n [rpm]'); xlabel(ax, 'timp [s]');

% 7. Eroarea efectorului față de r(t)
if isfield(res, 'eroare_r')
    ax = nexttile(tab_grafice(et, 'Eroare efector', 1, 1));
    plot(ax, t, 1000*res.eroare_r, 'LineWidth', 1.5); grid(ax, 'on');
    xlabel(ax, 'timp [s]'); ylabel(ax, '|r - r_d| [mm]');
    title(ax, 'Eroarea de poziție a efectorului față de traiectoria r(t)');
end
end

function legenda_tl(ax, etichete, titlu)
% O singură legendă, în dreapta grilei de grafice (nu acoperă curbele)
lg = legend(ax, etichete);
lg.Layout.Tile = 'east';
if nargin >= 3, title(lg, titlu); end
end

function linii3(ax, t, Y, cul, lw)
% Câte o linie pentru fiecare cuplă (culori fixe pentru 1, 2, 3)
for k = 1:3
    plot(ax, t, Y(k,:), 'Color', cul(k,:), 'LineWidth', lw);
end
end

%% ========================= VIZUALIZARE =================================
function tl = tab_grafice(etapa, nume, m, n)
% Toate graficele stau într-o singură fereastră "PUMA - Grafice": câte un tab
% pentru fiecare etapă și, în fiecare etapă, câte un sub-tab pentru fiecare grafic.
% Întoarce o grilă de m x n grafice (tiledlayout) în sub-tabul "nume".
fig = findobj(0, 'Type', 'figure', 'Tag', 'PUMA_grafice');
if isempty(fig)
    fig = figure('Name', 'PUMA - Grafice', 'NumberTitle', 'off', 'Tag', 'PUMA_grafice', ...
                 'Color', 'w', 'Position', [120 60 1150 720]);
    uitabgroup(fig, 'Tag', 'etape');
end
grup = findobj(fig(1), 'Tag', 'etape');
te = findobj(grup, '-depth', 1, 'Type', 'uitab', 'Title', etapa);
if isempty(te)
    te = uitab(grup, 'Title', etapa);
    uitabgroup(te, 'Tag', 'grafice');
    grup.SelectedTab = te;                        % se vede etapa curentă
end
sub = findobj(te(1), 'Tag', 'grafice');
delete(findobj(sub, '-depth', 1, 'Type', 'uitab', 'Title', nume));   % fără dubluri
tl = tiledlayout(uitab(sub, 'Title', nume), m, n, 'TileSpacing', 'compact', 'Padding', 'compact');
end

function h = vizualizare_robot(q, par, Rdorit, nume)
% Fereastra "Braț Robotic - Vizualizare Complexă": una singură, cu câte un tab
% pentru fiecare etapă (poziția statică și fiecare animație). În fiecare tab:
%   stânga  - vedere 3D
%   dreapta - proiecție în planul XY (de sus) și în planul YZ (din lateral)
% Robotul este desenat din cilindri (elemente albastre, cuple roșii),
% cu sistemul de coordonate al efectorului {e}. Unitățile afișate sunt cm.
if nargin < 3, Rdorit = []; end
if nargin < 4 || isempty(nume), nume = 'Robot'; end
fig = findobj(0, 'Type', 'figure', 'Tag', 'PUMA_robot');
if isempty(fig)
    fig = figure('Name', 'Braț Robotic - Vizualizare Complexă', 'NumberTitle', 'off', ...
                 'Tag', 'PUMA_robot', 'Color', 'w', 'Position', [80 80 1250 720]);
    uitabgroup(fig, 'Tag', 'animatii', 'SelectionChangedFcn', @(src, ~) opreste_animatii(src));
end
fig = fig(1);
grup = findobj(fig, 'Tag', 'animatii');
opreste_animatii(grup);
delete(findobj(grup, '-depth', 1, 'Type', 'uitab', 'Title', nume));   % fără dubluri
tab = uitab(grup, 'Title', nume, 'BackgroundColor', 'w');
grup.SelectedTab = tab;
figure(fig);
h.fig = fig;  h.tab = tab;
pozitii = {[0.03 0.12 0.52 0.76], [0.61 0.55 0.37 0.33], [0.61 0.12 0.37 0.33]};
for k = 1:3
    h.ax(k) = axes('Parent', tab, 'Position', pozitii{k});
end
titluri = {'Vedere 3D', 'Proiecție Plan XY (De sus)', 'Proiecție Plan YZ (Din lateral)'};
vederi  = {[40 25], [0 90], [90 0]};

piese = piese_robot(par);
L = 100*(par.l3 + par.l4 + par.le) + 3;                       % raza zonei de lucru [cm]
zmin = min(0, 100*(par.l1 + par.l2 - par.l4 - par.le)) - 3;
zmax = 100*(par.l1 + par.l2 + par.l4 + par.le) + 3;

for k = 1:3
    ax = h.ax(k);
    hold(ax, 'on');
    h.obiecte{k} = creeaza_piese(ax, piese);
    h.urma(k) = plot3(ax, NaN, NaN, NaN, 'm', 'LineWidth', 2);
    if ~isempty(Rdorit)
        plot3(ax, 100*Rdorit(1,:), 100*Rdorit(2,:), 100*Rdorit(3,:), 'k:', 'LineWidth', 1.2);
    end
    axis(ax, 'equal');
    axis(ax, [-L L -L L zmin zmax]);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'Axa X'); ylabel(ax, 'Axa Y'); zlabel(ax, 'Axa Z');
    title(ax, titluri{k});
    view(ax, vederi{k});
    light('Parent', ax, 'Position', [1 -1 2], 'Style', 'infinite');
    light('Parent', ax, 'Position', [-1 1 1], 'Style', 'infinite');
end
h.piese = piese;
h.titlu = uicontrol(tab, 'Style', 'text', 'Units', 'normalized', 'Position', [0 0.935 1 0.05], ...
                    'BackgroundColor', 'w', 'FontSize', 13, 'FontWeight', 'bold', 'String', '');
actualizare_robot(h, q, par);
tab.UserData = h;
end

function opreste_animatii(grup)
% Oprește orice animație care rulează în fereastra robotului
for tb = reshape(findobj(grup, '-depth', 1, 'Type', 'uitab'), 1, [])
    if isappdata(tb, 'ruleaza'), setappdata(tb, 'ruleaza', false); end
end
end

function actualizare_robot(h, q, par, P, tt)
% Mută robotul în configurația q (și desenează urma efectorului P, dacă există)
[~, ~, T] = geometric(q, par, 0);
cm = @(Tm) [Tm(1:3,1:3), 100*Tm(1:3,4); 0 0 0 1];
TT = {eye(4), cm(T.T01), cm(T.T02), cm(T.T03), cm(T.T0e)};   % cadrele 0, 1, 2, 3, e
for k = 1:numel(h.ax)
    ob = h.obiecte{k};
    for j = 1:numel(h.piese)
        pc = h.piese(j);
        Tm = TT{pc.cadru + 1};
        pts = Tm(1:3,1:3)*[pc.X(:).'; pc.Y(:).'; pc.Z(:).'] + Tm(1:3,4);
        sz = size(pc.X);
        if strcmp(pc.tip, 'text')
            set(ob(j), 'Position', pts.');
        else
            set(ob(j), 'XData', reshape(pts(1,:), sz), 'YData', reshape(pts(2,:), sz), ...
                       'ZData', reshape(pts(3,:), sz));
        end
    end
    if nargin >= 4 && ~isempty(P)
        set(h.urma(k), 'XData', 100*P(1,:), 'YData', 100*P(2,:), 'ZData', 100*P(3,:));
    end
end
qg = q*180/pi;
txt = sprintf('Braț Robotic: q1=%.1f° | q2=%.1f° | q3=%.1f°', qg(1), qg(2), qg(3));
if nargin >= 5 && ~isempty(tt)
    txt = sprintf('%s | t = %.2f s', txt, tt);
end
set(h.titlu, 'String', txt);
end

function piese = piese_robot(par)
% Geometria robotului (în cm), fiecare piesă descrisă în cadrul ei:
%   cadru 0 = baza fixă, 1 = {1}, 2 = {2}, 3 = {3}, 4 = {e}
s = 100;
r  = 100*par.raza;                      % raza elementelor [cm]
rc = 1.35*r;                            % raza cuplelor
albastru = [0.15 0.25 0.85];
rosu     = [0.85 0.15 0.15];
l1 = s*par.l1; l2 = s*par.l2; l3 = s*par.l3; l4 = s*par.l4; le = s*par.le;

piese = struct('cadru', {}, 'tip', {}, 'X', {}, 'Y', {}, 'Z', {}, 'culoare', {}, 'text', {});
piese = adauga_cilindru(piese, 0, [0 0 0], [0 0 0.6], 1.8*r, rosu);        % placa de bază
piese = adauga_cilindru(piese, 0, [0 0 0], [0 0 l1-rc], r, albastru);      % coloana bazei
piese = adauga_cilindru(piese, 1, [0 0 -rc], [0 0 rc], rc, rosu);          % cupla 1 (axa z)
piese = adauga_cilindru(piese, 1, [0 0 rc], [0 0 l2], r, albastru);        % element 1 - vertical
piese = adauga_sfera(piese, 1, [0 0 l2], r, albastru);                      % colț element 1
piese = adauga_cilindru(piese, 1, [0 0 l2], [l3 0 l2], r, albastru);       % element 1 - orizontal
piese = adauga_cilindru(piese, 2, [0 -rc 0], [0 rc 0], rc, rosu);          % cupla 2 (axa y)
piese = adauga_cilindru(piese, 2, [0 0 0], [l4 0 0], r, albastru);         % element 2
piese = adauga_cilindru(piese, 3, [0 -rc 0], [0 rc 0], rc, rosu);          % cupla 3 (axa y)
piese = adauga_cilindru(piese, 3, [0 0 0], [le 0 0], 0.8*r, albastru);     % element 3

% sistemul de coordonate al efectorului {e}
La = max(0.5*le, 4);
axe = {[La 0 0], [0 La 0], [0 0 La]};
culori = {[0.9 0 0], [0 0.7 0], [0 0 0.9]};
etichete = {'X', 'Y', 'Z'};
for k = 1:3
    v = axe{k};
    piese(end+1) = struct('cadru', 4, 'tip', 'linie', 'X', [0 v(1)], 'Y', [0 v(2)], ...
                          'Z', [0 v(3)], 'culoare', culori{k}, 'text', '');
    piese(end+1) = struct('cadru', 4, 'tip', 'text', 'X', 1.15*v(1), 'Y', 1.15*v(2), ...
                          'Z', 1.15*v(3), 'culoare', culori{k}, 'text', etichete{k});
end
end

function piese = adauga_cilindru(piese, cadru, p1, p2, r, culoare)
% Cilindru (cu capace) între p1 și p2, în cadrul "cadru"
n = 24;
th = linspace(0, 2*pi, n+1);
d = p2(:) - p1(:);  Lc = norm(d);  w = d/Lc;
if abs(w(1)) < 0.9, a = [1; 0; 0]; else, a = [0; 1; 0]; end
u = cross(a, w); u = u/norm(u);  v = cross(w, u);
Rm = [u v w];
% manta + capace: rândurile = (centru jos, margine jos, margine sus, centru sus)
rr = [0; r; r; 0] * ones(1, n+1);
zz = [0; 0; Lc; Lc] * ones(1, n+1);
xl = rr .* (ones(4,1)*cos(th));
yl = rr .* (ones(4,1)*sin(th));
pts = Rm*[xl(:).'; yl(:).'; zz(:).'] + p1(:);


piese(end+1) = struct('cadru', cadru, 'tip', 'suprafata', 'X', reshape(pts(1,:), 4, n+1), ...
                      'Y', reshape(pts(2,:), 4, n+1), 'Z', reshape(pts(3,:), 4, n+1), ...
                      'culoare', culoare, 'text', '');
end

function piese = adauga_sfera(piese, cadru, c, r, culoare)
[X, Y, Z] = sphere(16);
piese(end+1) = struct('cadru', cadru, 'tip', 'suprafata', 'X', r*X + c(1), 'Y', r*Y + c(2), ...
                      'Z', r*Z + c(3), 'culoare', culoare, 'text', '');
end

function ob = creeaza_piese(ax, piese)
ob = gobjects(1, numel(piese));
for j = 1:numel(piese)
    pc = piese(j);
    switch pc.tip
        case 'suprafata'
            ob(j) = surface(ax, pc.X, pc.Y, pc.Z, 'FaceColor', pc.culoare, 'EdgeColor', 'none', ...
                            'FaceLighting', 'gouraud', 'AmbientStrength', 0.5);
        case 'linie'
            ob(j) = plot3(ax, pc.X, pc.Y, pc.Z, 'Color', pc.culoare, 'LineWidth', 2.5);
        case 'text'
            ob(j) = text(ax, pc.X, pc.Y, pc.Z, pc.text, 'Color', pc.culoare, ...
                         'FontWeight', 'bold', 'FontSize', 11);
    end
end
end

function animatie(Q, par, Rdorit, viteza, t, nume)
% Animația robotului, într-un tab al ferestrei de vizualizare, cu butoane:
%   Redă (de unde a rămas / de la început), Pauză, bara de timp, viteza.
%   Q - 3 x N coordonatele cuplelor; Rdorit - traiectoria dorită (opțional)
%   viteza - 1 = timp real; t - timpul [s]; nume - titlul tabului
% Animația se redă o dată automat; apoi se poate relua oricând din tab.
if nargin < 3, Rdorit = []; end
if nargin < 4 || isempty(viteza), viteza = 1; end
N = size(Q, 2);
if nargin < 5 || isempty(t), t = 0.02*(0:N-1); end
if nargin < 6, nume = 'Animație'; end
P = zeros(3, N);
for i = 1:N
    P(:,i) = geometric(Q(:,i), par, 0);
end
h = vizualizare_robot(Q(:,1), par, Rdorit, nume);
tab = h.tab;
h.Q = Q;  h.t = t(:).';  h.P = P;  h.par = par;
h.viteze = [0.25 0.5 1 2 4];
[~, iv] = min(abs(h.viteze - viteza));
uicontrol(tab, 'Style', 'pushbutton', 'String', 'Redă', 'FontSize', 11, 'Units', 'normalized', ...
          'Position', [0.03 0.015 0.07 0.05], 'Callback', @(~, ~) reda(tab));
uicontrol(tab, 'Style', 'pushbutton', 'String', 'Pauză', 'FontSize', 11, 'Units', 'normalized', ...
          'Position', [0.11 0.015 0.07 0.05], 'Callback', @(~, ~) setappdata(tab, 'ruleaza', false));
h.bara = uicontrol(tab, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.20 0.025 0.56 0.035], ...
                   'Min', 1, 'Max', max(N, 2), 'Value', 1, 'SliderStep', [1 10]/max(N - 1, 1), ...
                   'Callback', @(s, ~) muta_bara(tab, s.Value));
try
    addlistener(h.bara, 'ContinuousValueChange', @(s, ~) muta_bara(tab, s.Value));
catch
end
uicontrol(tab, 'Style', 'text', 'String', 'viteza:', 'BackgroundColor', 'w', 'Units', 'normalized', ...
          'Position', [0.77 0.015 0.05 0.04]);
h.viteza = uicontrol(tab, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.82 0.02 0.07 0.04], ...
                     'String', {'0.25x', '0.5x', '1x', '2x', '4x'}, 'Value', iv);
h.timp = uicontrol(tab, 'Style', 'text', 'BackgroundColor', 'w', 'FontSize', 11, 'Units', 'normalized', ...
                   'Position', [0.90 0.015 0.09 0.04], 'String', sprintf('t = %.2f s', h.t(1)));
tab.UserData = h;
setappdata(tab, 'k', 1);
setappdata(tab, 'ruleaza', false);
reda(tab);
end

function reda(tab)
% Redă animația din tab în timp real (de la început, dacă s-a terminat)
if ~isvalid(tab) || getappdata(tab, 'ruleaza'), return; end
h = tab.UserData;
N = numel(h.t);
k = getappdata(tab, 'k');
if k >= N, k = 1; end
setappdata(tab, 'ruleaza', true);
tc = h.t(k);  ceas = tic;  ultim = 0;
while isvalid(tab) && getappdata(tab, 'ruleaza')
    acum = toc(ceas);
    tc = tc + (acum - ultim)*h.viteze(h.viteza.Value);  ultim = acum;
    k = min(N, max(k, find(h.t <= tc, 1, 'last')));
    arata_cadru(tab, k);
    if k >= N, break; end
    pause(0.005);
end
if isvalid(tab), setappdata(tab, 'ruleaza', false); drawnow; end
end

function muta_bara(tab, valoare)
% Bara de timp: oprește redarea și arată momentul ales
if ~isvalid(tab), return; end
setappdata(tab, 'ruleaza', false);
arata_cadru(tab, round(valoare));
drawnow;
end

function arata_cadru(tab, k)
% Desenează robotul în cadrul k al animației din tab
if ~isvalid(tab), return; end
h = tab.UserData;
k = max(1, min(numel(h.t), k));
setappdata(tab, 'k', k);
actualizare_robot(h, h.Q(:,k), h.par, h.P(:,1:k), h.t(k));
h.bara.Value = k;
h.timp.String = sprintf('t = %.2f s', h.t(k));
drawnow limitrate;
end
