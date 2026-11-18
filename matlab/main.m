clc
clear
close all

%% System Parameters and Maxrixes

l  = 0.235; % m
Ix = 8.11858 * 10^(-5); % kg/m^2
Iy = 8.11858 * 10^(-5); % kg/m^2
Iz = 6.12233 * 10^(-5); % kg/m^2
b1 = l/Ix;
b2 = l/Iy;
b3 = l/Iz;
m  = 0.506; % kg
g  = 9.8;   % m/s^2

A = [0 1 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 1 0 0;
     0 0 g 0 0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0 0 0 0 1;
    -g 0 0 0 0 0 0 0 0 0 0 0];

B = [0 0 0 0;
     0 b1 0 0; 
     0 0 0 0; 
     0 0 b2 0; 
     0 0 0 0; 
     0 0 0 b3; 
     0 0 0 0; 
     (1/m) 0 0 0; 
     0 0 0 0; 
     0 0 0 0;
     0 0 0 0;
     0 0 0 0];

C = [1 0 0 0 0 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0 0 0 0;
     0 0 0 0 0 0 1 0 0 0 0 0;
     0 0 0 0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 0 0 0 1 0];

D = 0; 

C1 = [0 0 0 0 1 0 0 0 0 0 0 0;
      0 0 0 0 0 0 1 0 0 0 0 0;
      0 0 0 0 0 0 0 0 1 0 0 0;
      0 0 0 0 0 0 0 0 0 0 1 0];
 
F = [[-1 1.33 0 0 0 0 0 0 0 0 0 0] ;
     [-1.33 -1 0 0 0 0 0 0 0 0 0 0];
     [0 0 -4 0 0 0 0 0 0 0 0 0];
     [0 0 0 -4 0 0 0 0 0 0 0 0];
     [0 0 0 0 -5 0 0 0 0 0 0 0];
     [0 0 0 0 0 -5 0 0 0 0 0 0];
     [0 0 0 0 0 0 -6 0 0 0 0 0];
     [0 0 0 0 0 0 0 -6 0 0 0 0];
     [0 0 0 0 0 0 0 0 -7 0 0 0];
     [0 0 0 0 0 0 0 0 0 -7 0 0];
     [0 0 0 0 0 0 0 0 0 0 -8 0];
     [0 0 0 0 0 0 0 0 0 0 0 -8]];
    
%% First part

sigma_v_n = 0.005;
sigma_w_d = 0.01;
R = eye(4);
v = [10, 10, 10, 10, 10, 10, 10, 10, 8, 8, 8, 8];
Q = diag(v);

syms s x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 t
syms phi phi_dot theta theta_dot psi psi_dot z z_dot x x_dot y y_dot

X = [phi phi_dot theta theta_dot psi psi_dot z z_dot x x_dot y y_dot];
a1 = (Iy - Iz) / Ix;
a2 = (Iz - Ix) / Iy;
a3 = (Ix - Iy) / Iz;
b1 = l / Ix;
b2 = l / Iy;
b3 = l / Iz;
u_t = cos(X(1))*cos(X(3));
u_x = cos(X(1))*sin(X(3))*cos(X(5)) + sin(X(1))*sin(X(5));
u_y = cos(X(1))*sin(X(3))*sin(X(5)) + sin(X(1))*cos(X(5));

sys = ss(A, B, C, D);

Gs = C * ((s * (eye(12)) - A) \ B) + D;

G = [0 tf([1.9079*10^(-5)],[1 0 0]) 0 0;
     0 0 tf([1.9079*10^(-5)],[1 0 0]) 0;
     0 0 0 tf([1.4387*10^(-5)],[1 0 0]);
     tf([1.9763],[1 0 0]) 0 0 0;
     0 0  tf([1.8697*10^(-4)],[1 0 0 0 0]) 0;
     0 tf([-1.8697*10^(-4)],[1 0 0 0 0]) 0 0;];
 
P = pole(G);
figure
pzmap(G)

R_A = rank(A);
R_B = rank(B);

disp('Controlability Matrix');
cm = ctrb(A,B);
disp('Rank of Controlability Matrix:');
rank(cm);

disp('rref');
rref(cm);

disp('Observability matrix')
obs = obsv(A,C)
disp('Rank of Observability Matrix:');
rank(obs);

% Compute the Jordan form
[V, J] = jordan(A);

% Display the Jordan form
disp("Jordan form:")
disp(J)

ej = [[1, t, t^2/2, t^3/6, 0, 0, 0, 0, 0, 0, 0, 0];
      [0, 1,   t,   t^2/2, 0, 0, 0, 0, 0, 0, 0, 0];
      [0, 0,   1,     t,   0, 0, 0, 0, 0, 0, 0, 0];
      [0, 0,   0,     1,   0, 0, 0, 0, 0, 0, 0, 0];
      [0, 0,   0,     0,   1, t, 0, 0, 0, 0, 0, 0];
      [0, 0,   0,     0,   0, 1, 0, 0, 0, 0, 0, 0];
      [0, 0,   0,     0,   0, 0, 1, t, 0, 0, 0, 0];
      [0, 0,   0,     0,   0, 0, 0, 1, 0, 0, 0, 0];
      [0, 0,   0,     0,   0, 0, 0, 0, 1, t, 0, 0];
      [0, 0,   0,     0,   0, 0, 0, 0, 0, 1, 0, 0];
      [0, 0,   0,     0,   0, 0, 0, 0, 0, 0, 1, t];
      [0, 0,   0,     0,   0, 0, 0, 0, 0, 0, 0, 1]];


% Display the similarity transformation matrix
disp("Similarity transformation matrix:")
disp(V)
rank(V)

V1 = inv(V);     % V inverse

ea = V * ej / V;     % V * e^j * inv(V)

x0 = [x1; x2; x3; x4; x5; x6; x7; x8; x9; x10; x11; x12];
 
c_Q7 = V \ x0;     % inv(V) * x0

T =  [0,       0, 1.0000,      0,      0,      0,      0,      0,      0,      0,       0,       0;
      0,       0,      0, 1.0000,      0,      0,      0,      0,      0,      0,       0,       0;
      0,       0,      0,      0,      0,      0,      0,      0, 1.0000,      0,       0,       0;
      0,       0,      0,      0,      0,      0,      0,      0,      0, 1.0000,       0,       0;
      0,       0,      0,      0, 1.0000,      0,      0,      0,      0,      0,       0,       0;
      0,       0,      0,      0,      0, 1.0000,      0,      0,      0,      0,       0,       0;
      0,       0,      0,      0,      0,      0, 1.0000,      0,      0,      0,       0,       0;
      0,       0,      0,      0,      0,      0,      0, 1.0000,      0,      0,       0,       0;
      0,  9.8000,      0,      0,      0,      0,      0,      0,      0,      0, -9.8000,       0;
      0,       0, 9.8000,      0,      0,      0,      0,      0,      0,      0,       0, -9.8000;
-9.8000,       0,      0,      0,      0,      0,      0,      0,      0,      0,       0,       0;
      0, -9.8000,      0,      0,      0,      0,      0,      0,      0,      0,       0,       0];


figure
step(sys)
x0 = [1;0;0;1;0;1;0;1;0;1;0;1];
initial(sys,x0)

%% Second part

xo = [1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

eig(F)
                  
k = place(A, B, [-6.8-7.2i, -6.8+7.2i, -35, -35, -36, -36, -37, -37, -38, -38, -39, -39]);

eigen_value = eig(A);

poles = [-6.8-7.2i,-6.8+7.2i,-35,-35,-36,-36,-37,-37,-38,-38,-39,-39];
% poles = [-1-1.33i, -1+1.33i, -4, -4, -5, -5, -6, -6, -7, -7, -8, -8];
K_place = place(A,B,poles);

% Static Pre-compensator
A_cl = A - B*K_place;
G_CL_0 = -C1 / (A_cl) * B;  % -C1 * inv(A_cl) * B
inv_G_CL_0 = inv(G_CL_0);

% Dynamic Pre-compensator
A_bar = [A zeros(12,4); -C1 zeros(4,4)];
B_bar = [B; zeros(4,4)];
Q_bar = diag([10, 5, 10, 5, 10, 10, 50, 10, 50, 100, 100, 100, 100, 100, 100, 100]);
Q_bar1 = diag([10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]);
R = 1;

[K_LQR, ~, ~] = lqr(A_bar, B_bar, Q_bar1, R, 0);
K_d_LQR = K_LQR(:, 1:12);
K_q_LQR = K_LQR(:, 13:16);

% syms k_1_1 k_1_2 k_1_3 k_1_4 k_1_5 k_1_6 k_1_7 k_1_8 k_1_9 k_1_10 k_1_11 k_1_12 k_2_1 k_2_2 k_2_3 k_2_4 k_2_5 k_2_6 k_2_7 k_2_8 k_2_9 k_2_10 k_2_11 k_2_12 k_3_1 k_3_2 k_3_3 k_3_4 k_3_5 k_3_6 k_3_7 k_3_8 k_3_9 k_3_10 k_3_11 k_3_12 k_4_1 k_4_2 k_4_3 k_4_4 k_4_5 k_4_6 k_4_7 k_4_8 k_4_9 k_4_10 k_4_11 k_4_12
% syms k_q11 k_q12 k_q13 k_q14 k_q21 k_q22 k_q23 k_q24 k_q31 k_q32 k_q33 k_q34 k_q41 k_q42 k_q43 k_q44 s
% 
% K_q = [k_q11 k_q12 k_q13 k_q14;
%        k_q21 k_q22 k_q23 k_q24;
%        k_q31 k_q32 k_q33 k_q34;
%        k_q41 k_q42 k_q43 k_q44];
% K_d = [k_1_1 k_1_2 k_1_3 k_1_4 k_1_5 k_1_6 k_1_7 k_1_8 k_1_9 k_1_10 k_1_11 k_1_12;
%        k_2_1 k_2_2 k_2_3 k_2_4 k_2_5 k_2_6 k_2_7 k_2_8 k_2_9 k_2_10 k_2_11 k_2_12;
%        k_3_1 k_3_2 k_3_3 k_3_4 k_3_5 k_3_6 k_3_7 k_3_8 k_3_9 k_3_10 k_3_11 k_3_12;
%        k_4_1 k_4_2 k_4_3 k_4_4 k_4_5 k_4_6 k_4_7 k_4_8 k_4_9 k_4_10 k_4_11 k_4_12];
% K_D_S = [K_d K_q];
% 
% if (rank([B A; zeros(4,4) -C1]) == length(Ai))  % Controllability Condition for Dynamic P.Comp.
%     disp ('Augmented system with Integral State is controllable')
%     A_CL_D = [A-(B*K_d) B*K_q;
%               -C1 zeros(4, 4)]
%     CE_D = vpa(det(s*eye(16)-A_CL_D), 16);
%     vpa(collect (CE_D,s),12)
%     coeff_D = coeffs(CE_D,s);
%     coeff_D = vpa(fliplr(coeff_D),12)
%     alpha_D = ((s+4)^2)*((s+5)^2)*((s+6)^2)*((s+7)^2)*((s+8)^2)*(s^2+2*s+2.7689); % Desired Poles
%     coeff_alpha = coeffs(alpha_D,s);
%     coeff_alpha = vpa(fliplr(coeff_alpha),12)
%     for i = 1:length(A_CL_D)
%         equns(i)=coeff_D(i+1)-coeff_alpha(i+1);
%     end
%     K_D = solve(equns,[k_1 k_2 k_3 k_4 k_5 k_6 k_7 k_8 k_9 k_10 k_11 k_12 K_q]);
%     K_d = double([K_D.k_1 K_D.k_2 K_D.k_3 K_D.k_4 K_D.k_5 K_D.k_6 K_D.k_7 K_D.k_8 K_D.k_9 K_D.k_10 K_D.k_11 K_D.k_12])
%     K_q = double(K_D.k_q)
% else
%     disp ('Augmented system with Integral State is not controllable')
% end
% 

%% third part

disp('observer (full rank)')
% AT-CT*K & L=KT
AT = transpose(A);
CT = transpose(C1);
% det(AT-CT*K)=0 should gives us the desired poles
% our desired poles are:
poles = [-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13];
%[C, K] = place(A, poles);
K = place(AT,CT, poles);
L = K';

% Reduced observer
disp('observer (reduced)')
n = 12;
q = rank(C1); % q=4 & n-q=8
F_r = diag([-2,-2,-3,-3,-4,-4,-5,-5]);
L_r = [rand(8,4)];
rank(ctrb(F_r, L_r)); % FullRank = 8
K_r = transpose(L_r);

%Calculating T through the solution of the lyapanov function
T = lyap(F_r, -A, L_r * C1)
inv_c_T = inv([C1; T]);
Q_1 = inv_c_T(:,1:q)
Q_2 = inv_c_T(:,q+1:n)
