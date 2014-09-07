% (C) 2012 PRL Dev Team. All rights reserved.
%   Author: Przemek Lekston <p.lekston@gmail.com>

% Simple Rocket control simulation for a 3 DoF model of a lateral motion of
% a vertical rocket propelled by a reaction engine mounted at its bottom. 
% (embedded in 2D space)
% Rocket state is defined by:

%  theta_dot    - tilt angular rate
%  theta        - tilt (zero when vertical)
%  x_dot        - lateral velocity

% The controller works best at stabilizing the lateral speed of the rocket (3rd variable)

clear all;

c = 0.01;   %N/(m/s) - a very minimalistic representation of air drag
m = 10;     %kg
len = 2;      %m
g = 9.81;    %m/s^2

%gains for the feedback controller
k_th = 500;
k_th_dot = 200;
k_x_dot = -80;

%20; 0.5; -0.5

t = 0:0.01:6;

%demand (require 
u = [0*ones(size(t));0*ones(size(t));8*ones(size(t))]';

k_th_dot_prim = k_th_dot + k_x_dot*len

s = tf('s');

Den = [1, -c/m, -g/len, len/m]; %coeffs of Den in descending powers of s

Num = [ 0, 1/(m*len), -(c/m)/(m*len) ; 
        0, 1/m*len, -1/m*len*c/m     ; 
        (1-len)/m,  -1/m*len*c*len/m, - 1/m*len*g*len - (g/len)*(1-len)/m];

[rock.a,rock.b,rock.c,rock.d] = tf2ss(Num,Den);

%alternatively:
%Num = { [0, 1/(m*len), -(c/m)/(m*len)] ; 
%        [0, 1/m*len, -1/m*len*c/m]     ; 
%        [(1-len)/m,  -1/m*len*c*len/m, - 1/m*len*g*len - (g/len)*(1-len)/m]};
%
%rock_tf = tf(Num, {Den ; Den ; Den});

rock_ss = ss(rock.a,rock.b,rock.c,rock.d);


Ctrl_Num = [{k_th_dot_prim}, {k_th}, {k_x_dot}];

CtrlTF = tf(Ctrl_Num,[{1},{1},{1}]);    %3 inputs, 1 output

SensorTF = tf(eye(3));     %sensor impact disregarded

C_loop = feedback(series(CtrlTF,rock_ss),SensorTF,[1,2,3],[1,2,3]);

lsim(C_loop,u,t);

[A,B,C,D] = ssdata(C_loop);

[N,D] = ss2tf(A,B,C,D,1);

roots(D)
