clear all
close all
clc

%% Part A

g = 9.81;
R = 2;
L = 10^(-6);
K = 10^(-5);
M = 0.020;
y0 = 0.01;
i0 = y0*sqrt((M*g)/K);
u0 = R*i0;
% Create state space matrix
A = [0 1 0;
    2*g/y0 0 -2*sqrt((K*g)/(M*y0^2));
    (1/L)*(u0-R*i0) 0 -R*y0/L];
B = [0;0;y0/L];
C = [1 0 0];
D = [0];
% Create state space model
Gc_ss = ss(A,B,C,D);
[b,a] = ss2tf(A,B,C,D);
Gc = tf(b,a);
% Convert to discrete
Ts = 0.001;
Gd_ss = c2d(Gc_ss,Ts);
G = Gd_ss.A;
H = Gd_ss.B;
Gd = tf(Gd_ss);
% (-1) gain for reverse KVL
Gd = -Gd;

%% Part B

kisi = 0.69;
wn = 18;
s1 = -kisi*wn + 1i*wn*sqrt(1-kisi^2);
s2 = -kisi*wn - 1i*wn*sqrt(1-kisi^2);
z1 = exp(s1*Ts);
z2 = exp(s2*Ts);
z = tf('z',Ts);
e_ss = 10;
y00 = 0.012;
Dd = 487.25*(z-0.9567)/(z-0.9315);
% stepinfo(feedback(Gd*Dd,1))
% step(feedback(Gd*Dd,1))

%% Part C

[T,eigV] = eig(G);
dG = T^(-1)*G*T;
dH = T^(-1)*H;
dC = C*T;
syms zz;
Mn = [H G*H G*G*H];
a_coeff = fliplr(double(coeffs(det(zz*eye(3)-G))));
W = [a_coeff(3),a_coeff(2),1;
    a_coeff(2),1,0;
    1,0,0];
T = Mn*W;
ddG = T^(-1)*G*T;
ddH = T^(-1)*H;

z3 = 0.5;
alpha_coeff = fliplr(double(coeffs((zz-z1)*(zz-z2)*(zz-z3))));
Kc = [(alpha_coeff(4)-a_coeff(4)),(alpha_coeff(3)-a_coeff(3)),(alpha_coeff(2)-a_coeff(2))]*T^(-1);

Phi_G = G^(3);
Lo = Phi_G*[C;C*G;C*G*G]^(-1)*[0;0;1];

N = (-C*(G-eye(3)-H*Kc)^(-1)*H)^(-1);

%% Part D

Ad = 0.02;
wd = 2*pi*10;

% Create state space matrix
A_d = [0 1;
    -wd^2 0];
B_d = [0;0];
C_d = [1 0];
D_d = [0];
% Create state space model
Gdc_ss = ss(A_d,B_d,C_d,D_d);
% Convert to discrete
Ts = 0.001;
Gdd_ss = c2d(Gdc_ss,Ts);
G_d = Gdd_ss.A;
H_d = Gdd_ss.B;
Gd = tf(Gdd_ss);

Gw = [G H*C_d;[0 0 0;0 0 0] G_d];
Hw = [H;[0;0]];
Cw = [C 0 0];
Dw = [0];

Phi_G_w = Gw^(5);
Low = Phi_G_w*[Cw;Cw*Gw;Cw*Gw*Gw;Cw*Gw*Gw*Gw;Cw*Gw*Gw*Gw*Gw]^(-1)*[0;0;0;0;1];
L_p = [Low(1);Low(2);Low(3)];
L_d = [Low(4);Low(5)];