clear all

%%%%%%%%%%%% Walce %%%%%%%%%%%%%
R1 = 1.12838;
H1 = 5;
c1 = 4.3589;

%%%%%%%%%%%% Ostroslup %%%%%%%%%%%%%
a2 = 3.59153;
b2 = 7.79612;
H2 = 19;
c2 = 3.08221;

%%%%%%%%%%%% Charakterystyka statyczna %%%%%%%%%%%%%
qwe = 0.5*H2;
h10 = (qwe/c1)^2;
h20 = (qwe/c2)^2;
h1i = 0;
h2i = 0.1*h20;
qdiff = sqrt(h20*1.1)*c2-qwe;

%%%%%%%%%%%% Parametry oraz transmitancje %%%%%%%%%%%%%
k = 2*sqrt(h20)/c2;
T1 = 2*pi*sqrt(h10)*R1^2/c1;
T2 = 2*a2*b2*h20^(5/2)/(H2^2*c2);
T3 = 0.2*min(T1, T2);

s = tf('s');
Kob=k/((1+T1*s)*(1+T2*s)*(1+T3*s));

%%
%%%%%%%%%%%% Regulator P %%%%%%%%%%%%%

figure
rlocus(Kob);
hold on
fplot(@(x)-2*x, [-5, 0]);
kr_P_gr = 14.3;
kr_P = 1.92; % wzmocnienie odczytane z wykresu
n_P = 0.128; % stopien stabilnosci odczytany z wykresu
K_P = kr_P*Kob;
G_P = K_P/(1+K_P);
figure
step(G_P, [0:0.001:70]);
grid on

%%
%%%%%%%%%%%% Regulator PI %%%%%%%%%%%%%

kr_PI_Basic = 1;
Tc_PI_x = T1;
Tc_PI_y = T2;
Tc_PI_z = T3;
K_PI_Basic_x = kr_PI_Basic*(1+Tc_PI_x*s)/s*Kob;
K_PI_Basic_y = kr_PI_Basic*(1+Tc_PI_y*s)/s*Kob;
K_PI_Basic_z = kr_PI_Basic*(1+Tc_PI_z*s)/s*Kob;
% figure
% rlocus(K_PI_Basic_x);
% rlocus(K_PI_Basic_y);
% rlocus(K_PI_Basic_z);
% hold on
% fplot(@(x)-2*x, [-5, 0]);
kr_PI_gr_x = margin(K_PI_Basic_x);
kr_PI_gr_y = margin(K_PI_Basic_y);
kr_PI_gr_z = margin(K_PI_Basic_z);
kr_PI_x = 0.0391;
kr_PI_y = 0.105;
kr_PI_z = 0.0267;
n_PI_x = 0.0334;
n_PI_y = 0.101;
n_PI_z = 0.0267;
K_PI_x = kr_PI_x*K_PI_Basic_x;
K_PI_y = kr_PI_y*K_PI_Basic_y;
K_PI_z = kr_PI_z*K_PI_Basic_z;
G_PI_x = minreal(K_PI_x/(1+K_PI_x));
G_PI_y = minreal(K_PI_y/(1+K_PI_y));
G_PI_z = minreal(K_PI_z/(1+K_PI_z));

% Parametry czasowe
stpinf_PI_x = stepinfo(G_PI_x);
stpinf_PI_y = stepinfo(G_PI_y);
stpinf_PI_z = stepinfo(G_PI_z);
eust_PI_x = 1 - dcgain(G_PI_x);
eust_PI_y = 1 - dcgain(G_PI_y);
eust_PI_z = 1 - dcgain(G_PI_z);
przeregulowanie_procent_PI_x = stpinf_PI_x.Overshoot;
przeregulowanie_procent_PI_y = stpinf_PI_y.Overshoot;
przeregulowanie_procent_PI_z = stpinf_PI_z.Overshoot;
czas_narastania_PI_x = stpinf_PI_x.RiseTime;
czas_narastania_PI_y = stpinf_PI_y.RiseTime;
czas_narastania_PI_z = stpinf_PI_z.RiseTime;
tr_PI_x = stepinfo(G_PI_x, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PI_y = stepinfo(G_PI_y, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PI_z = stepinfo(G_PI_z, 'SettlingTimeTreshold', 0.02).TransientTime;

figure
step(G_PI_x, G_PI_y, G_PI_z,  [0:0.01:250]);
legend("Tc=T1", "Tc=T2", "Tc=T3");
grid on

%%
%%%%%%%%%%%% Regulator PD %%%%%%%%%%%%%

kr_PD_Basic = 1;
Tb_PD = 0.1*min([T1, T2, T3]);
Ta_PD_x = T1;
Ta_PD_y = T2;
Ta_PD_z = T3;
K_PD_Basic_x = kr_PD_Basic*(1+Ta_PD_x*s)/(1+Tb_PD*s)*Kob;
K_PD_Basic_y = kr_PD_Basic*(1+Ta_PD_y*s)/(1+Tb_PD*s)*Kob;
K_PD_Basic_z = kr_PD_Basic*(1+Ta_PD_z*s)/(1+Tb_PD*s)*Kob;
% figure
% rlocus(K_PD_Basic_x);
% rlocus(K_PD_Basic_y);
% rlocus(K_PD_Basic_z);
% hold on
% fplot(@(x)-2*x, [-5, 0]);
kr_PD_gr_x = margin(K_PD_Basic_x);
kr_PD_gr_y = margin(K_PD_Basic_y);
kr_PD_gr_z = margin(K_PD_Basic_z);
kr_PD_x = 9.03;
kr_PD_y = 3.08;
kr_PD_z = 2.93;
n_PD_x = 0.588;
n_PD_y = 0.665;
n_PD_z = 0.156;
K_PD_x = kr_PD_x*K_PD_Basic_x;
K_PD_y = kr_PD_y*K_PD_Basic_y;
K_PD_z = kr_PD_z*K_PD_Basic_z;
G_PD_x = K_PD_x/(1+K_PD_x);
G_PD_y = K_PD_y/(1+K_PD_y);
G_PD_z = K_PD_z/(1+K_PD_z);

% Parametry czasowe
stpinf_PD_x = stepinfo(G_PD_x);
stpinf_PD_y = stepinfo(G_PD_y);
stpinf_PD_z = stepinfo(G_PD_z);
eust_PD_x = 1 - dcgain(G_PD_x);
eust_PD_y = 1 - dcgain(G_PD_y);
eust_PD_z = 1 - dcgain(G_PD_z);
przeregulowanie_procent_PD_x = stpinf_PD_x.Overshoot;
przeregulowanie_procent_PD_y = stpinf_PD_y.Overshoot;
przeregulowanie_procent_PD_z = stpinf_PD_z.Overshoot;
czas_narastania_PD_x = stpinf_PD_x.RiseTime;
czas_narastania_PD_y = stpinf_PD_y.RiseTime;
czas_narastania_PD_z = stpinf_PD_z.RiseTime;
tr_PD_x = stepinfo(G_PD_x, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PD_y = stepinfo(G_PD_y, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PD_z = stepinfo(G_PD_z, 'SettlingTimeTreshold', 0.02).TransientTime;

figure
step(G_PD_x, G_PD_y, G_PD_z, [0:0.001:70]);
legend("Tc=T1", "Tc=T2", "Tc=T3");
grid on

%%
%%%%%%%%%%%% Regulator PID %%%%%%%%%%%%%

kr_PID_Basic = 1;
Tb_PID = 0.1*min([T1, T2, T3]);
Ta_PID_x = T1;
Tc_PID_x = T2;
Ta_PID_y = T1;
Tc_PID_y = T3;
Ta_PID_z = T2;
Tc_PID_z = T3;
K_PID_Basic_x = minreal(kr_PID_Basic*(1+Ta_PID_x*s)*(1+Tc_PID_x*s)/(s*(1+Tb_PID*s))*Kob);
K_PID_Basic_y = minreal(kr_PID_Basic*(1+Ta_PID_y*s)*(1+Tc_PID_y*s)/(s*(1+Tb_PID*s))*Kob);
K_PID_Basic_z = minreal(kr_PID_Basic*(1+Ta_PID_z*s)*(1+Tc_PID_z*s)/(s*(1+Tb_PID*s))*Kob);
% figure
% rlocus(K_PID_Basic_x);
% rlocus(K_PID_Basic_y);
% rlocus(K_PID_Basic_z);
% hold on
% fplot(@(x)-2*x, [-5, 0]);
kr_PID_gr_x = margin(K_PID_Basic_x);
kr_PID_gr_y = margin(K_PID_Basic_y);
kr_PID_gr_z = margin(K_PID_Basic_z);
kr_PID_x = 0.627;
kr_PID_y = 0.044;
kr_PID_z = 0.149;
n_PID_x = 0.557;
n_PID_y = 0.0355;
n_PID_Z = 0.122;
K_PID_x = kr_PID_x*K_PID_Basic_x;
K_PID_y = kr_PID_y*K_PID_Basic_y;
K_PID_z = kr_PID_z*K_PID_Basic_z;
G_PID_x = K_PID_x/(1+K_PID_x);
G_PID_y = K_PID_y/(1+K_PID_y);
G_PID_z = K_PID_z/(1+K_PID_z);

% Parametry czasowe
stpinf_PID_x = stepinfo(G_PID_x);
stpinf_PID_y = stepinfo(G_PID_y);
stpinf_PID_z = stepinfo(G_PID_z);
eust_PID_x = 1 - dcgain(G_PID_x);
eust_PID_y = 1 - dcgain(G_PID_y);
eust_PID_z = 1 - dcgain(G_PID_z);
przeregulowanie_procent_PID_x = stpinf_PID_x.Overshoot;
przeregulowanie_procent_PID_y = stpinf_PID_y.Overshoot;
przeregulowanie_procent_PID_z = stpinf_PID_z.Overshoot;
czas_narastania_PID_x = stpinf_PID_x.RiseTime;
czas_narastania_PID_y = stpinf_PID_y.RiseTime;
czas_narastania_PID_z = stpinf_PID_z.RiseTime;
tr_PID_x = stepinfo(G_PID_x, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PID_y = stepinfo(G_PID_y, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PID_z = stepinfo(G_PID_z, 'SettlingTimeTreshold', 0.02).TransientTime;

figure
step(G_PID_x, G_PID_y, G_PID_z, [0:0.001:200]);
legend("Ta=T1, Tc=T2", "Ta=T1, Tc=T3", "Ta=T2, Tc=T3");
grid on

%%
%%%%%%%%%%%% Regulator PID niedokładny %%%%%%%%%%%%%

kr_PIDN_Basic = 1;
Tb_PIDN = 0.1*min([T1, T2, T3]);
Ta_PIDN_x = T1+0.1*T1;
Tc_PIDN_x = T2+0.1*T2;
Ta_PIDN_y = T1-0.1*T1;
Tc_PIDN_y = T2-0.1*T2;
K_PIDN_Basic_x = minreal(kr_PIDN_Basic*(1+Ta_PIDN_x*s)*(1+Tc_PIDN_x*s)/(s*(1+Tb_PIDN*s))*Kob);
K_PIDN_Basic_y = minreal(kr_PIDN_Basic*(1+Ta_PIDN_y*s)*(1+Tc_PIDN_y*s)/(s*(1+Tb_PIDN*s))*Kob);
figure
rlocus(K_PIDN_Basic_x);
rlocus(K_PIDN_Basic_y);
hold on
fplot(@(x)-2*x, [-5, 0]);
kr_PIDN_gr_x = margin(K_PIDN_Basic_x);
kr_PIDN_gr_y = margin(K_PIDN_Basic_y);
kr_PIDN_x = 0.534;
kr_PIDN_y = 0.744;
n_PIDN_x = 0.571;
n_PIDN_y = 0.538;
K_PIDN_x = kr_PIDN_x*K_PIDN_Basic_x;
K_PIDN_y = kr_PIDN_y*K_PIDN_Basic_y;
G_PIDN_x = K_PIDN_x/(1+K_PIDN_x);
G_PIDN_y = K_PIDN_y/(1+K_PIDN_y);

% Parametry czasowe
stpinf_PIDN_x = stepinfo(G_PIDN_x);
stpinf_PIDN_y = stepinfo(G_PIDN_y);
eust_PIDN_x = 1 - dcgain(G_PIDN_x);
eust_PIDN_y = 1 - dcgain(G_PIDN_y);
przeregulowanie_procent_PIDN_x = stpinf_PIDN_x.Overshoot;
przeregulowanie_procent_PIDN_y = stpinf_PIDN_y.Overshoot;
czas_narastania_PIDN_x = stpinf_PIDN_x.RiseTime;
czas_narastania_PIDN_y = stpinf_PIDN_y.RiseTime;
tr_PIDN_x = stepinfo(G_PIDN_x, 'SettlingTimeTreshold', 0.02).TransientTime;
tr_PIDN_y = stepinfo(G_PIDN_y, 'SettlingTimeTreshold', 0.02).TransientTime;

figure
step(G_PIDN_y, G_PID_x, G_PIDN_x, [0:0.001:20]);
legend("T-10%", "T", "T+10%");
grid on

















