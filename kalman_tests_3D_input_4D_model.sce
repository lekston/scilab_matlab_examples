// (C) 2012 PRL Dev Team. All rights reserved.
//   Author: Przemek Lekston <p.lekston@gmail.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in
//    the documentation and/or other materials provided with the
//    distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
// OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
// AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Kalman filter used for approximating the Aircraft sideslip angle 
// by using ban angle & yaw- & roll-rates with a given linear model of the lateral motion

clear;

//Generating Noise(semi-Gaussian)
T_0 = 10; //0.1Hz
T_1 = 2; //0.5Hz
T_2 = 0.2; //5Hz
T_3 = 0.0125; //80Hz

k_T0 = 0.02;
k_T1 = 0.03;
k_T2 = 0.05;
k_T3 = 0.4;

f_dscr = 50; //Hz

t = 0:(1/f_dscr):10;
v2 = zeros(1,size(t,2));

for i = 1:size(t,2)
v2(i) = sin(t(i)*2*%pi/T_0)*k_T0;
v2(i) = v2(i) + sin(t(i)*2*%pi/T_1)*k_T1;
v2(i) = v2(i) + sin(t(i)*2*%pi/T_2)*k_T2;
v2(i) = v2(i) + sin(t(i)*2*%pi/T_3)*k_T3;
end


// ******************************** A/C model **********************************
//operational conditions for this model:
//Altitude 0 ft
//Speed 30kts
A = [ -0.191397946662194	-0.777349563278821	-0.107877354174717	0.635321413922648;
       0.215155508053696	-0.373814775095028	-0.802673098904095	0;
      -0.823062092462765	 1.27100781273508	 -10.4944285156701	 0;
       0	                -0.0315276391144913    1                	0];
       
//             \beta         r         p      \phi
//   \beta   -0.1467   -0.7771  -0.03824    0.8285
//   r        0.1108   -0.3306    -1.033         0
//   p       -0.5262     2.776    -7.969         0
//   \phi          0   0.03816         1         0

B = [-0.0759529807611695	 0.104944810392342;
      0.0585428627522038	-0.132941251704353;
     -1.48714069050123	   0.00171011187011539;
      0	0];

//               \xi     \zeta
//   \beta  -0.06257   0.08181
//   r       0.03702  -0.07945
//   p       -0.9039  0.001022
//   \phi          0         0
 
C = [57.3      0      0      0;
      0       57.3    0      0;
      0        0     57.3    0;
      0        0      0   57.3];
//                \beta      r      p   \phi
//   \beta [deg]   57.3      0      0      0
//   r [deg/s]        0   57.3      0      0
//   p [deg/s]        0      0   57.3      0
//   \phi [deg]       0      0      0   57.3
 
D = [0 0;
     0 0;
     0 0;
     0 0];
//                  \xi  \zeta
//   \beta [deg]      0      0
//   r [deg/s]        0      0
//   p [deg/s]        0      0
//   \phi [deg]       0      0

X0 = [0*(%pi/180);
      1.7*(%pi/180);
      0*(%pi/180);
      3*(%pi/180)];
//quasi static at xi = 1deg, dzeta = -4.45

//Initial control deflections
xi0 = 6 * (%pi/180)
dzeta0 = 10 * (%pi/180)
          
SI = syslin('c', A, B, C, D, X0);

u = zeros(2, size(t,2));
for i = 1:size(t,2)
    u(:,i) = [xi0; dzeta0]; 
end

resp = csim(u,t,SI); 


// ************************** KALMAN FILTER ALGORITHM ********************

//initialization
SI_dscr = dscr(SI,(1/f_dscr));

A_dscr = SI_dscr(2);
B_dscr = SI_dscr(3);

X_est = X0;
X_KF = X0;
Inn = zeros(3,1);
S = zeros(3,3);
z = zeros(3,1);

//generating measurement noise
v3 = grand(1,size(t,2),'unf',-1,1);
v4 = grand(1,size(t,2),'unf',-1,1);

v(1,1,:) = v2(1,:); 
v(1,2,:) = zeros(1,size(t,2));
v(1,3,:) = zeros(1,size(t,2));

v(2,1,:) = zeros(1,size(t,2));
v(2,2,:) = v3(1,:);
v(2,3,:) = zeros(1,size(t,2));

v(3,1,:) = zeros(1,size(t,2));
v(3,2,:) = zeros(1,size(t,2));
v(3,3,:) = v3(1,:);

//computing measurement noise covariance
[R2, MeasMean] = corr(v2,1);
[R23, MeasMean] = corr(v2,v3,1);
[R24, MeasMean] = corr(v2,v4,1);

R32 = R23;
[R3, MeasMean] = corr(v3,1);
[R34, MeasMean] = corr(v3,v4,1);

R42 = R24;
R43 = R34;
[R4, MeasMean] = corr(v4,1);

//REVISION OF COVARIANCE MATRIX IS CRUCIAL HERE!!

R = [R2  R23 R24; 
     R32 R3  R34;
     R42 R43 R4;]

//Noise model (just to ensure proper alignemnt of matrices)
D = [1; 1; 1];

//assuming model noise covariance
Qii = R2/10 + R3/5 + R4/4;

//REVISION OF MODEL NOISE COVARIANCE MATRIX IS CRUCIAL HERE!!
Q = eye(4,4) * Qii;

//assuming initial value of the process noise
Pr_cov = eye(4,4) * 0.01;
Pr_cov_est = Pr_cov;

// since measurement is limited to bank angle, pitch rate and yaw rate; 
// the measurement matrix H != Id
H = [0 1 0 0;
     0 0 1 0;
     0 0 0 1];

//calculating filter output
for i = 1:(size(t,2)-1)
    
    //predict the state
    X_est(:,i+1) = A_dscr*X_KF(:,i) + B_dscr*[xi0; dzeta0]; 
    
    //predict the state covariance
    Pr_cov_est(:,:,i+1) = A_dscr*Pr_cov(:,:,i)*A_dscr' + Q; //4x4



    //noisy measurement
    z(:,i+1) = H*resp(:,i+1) + v(:,:,i+1)*D; //3x1

    
    //Innovation
    Inn(:,i+1) = z(:,i+1) - H*X_est(:,i+1); //3x1
    
    //Covariance of the innovation
    S(:,:,i+1) = R + H*Pr_cov_est(:,:,i+1)*H'; //3x3
    
    //update Kalman gain
    K(:,:,i+1) = Pr_cov_est(:,:,i+1)*H'*inv(S(:,:,i+1)); //4x3



    //update state
    X_KF(:,i+1) = X_est(:,i+1) + K(:,:,i+1)*Inn(:,i+1); //4x1
    
    //update the state covariance
    Pr_cov(:,:,i+1) = Pr_cov_est(:,:,i+1) - K(:,:,i+1)*S(:,:,i+1)*(K(:,:,i+1)'); //4x4
    
    
end


// ********************** END OF KALMAN FILTER ALGORITHM *****************

//figure 0 - sideslip beta
scf(0);
clf();
plot(t,resp(1,:),'b-');
plot(t,X_KF(1,:),'y-');
xtitle("sideslip");

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;

//figure 1 - yaw rate - r
scf(1);
clf();
plot(t,resp(2,:),'b-');
plot(t,z(1,:),'r-');
plot(t,X_KF(2,:),'y-');
xtitle("yaw rate");

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;

//figure 2 - roll rate - p
scf(2);
clf();
plot(t,resp(3,:),'b-');
plot(t,z(2,:),'r-');
plot(t,X_KF(3,:),'y-');
xtitle("roll rate");

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;

//figure 3 - bank angle - phi
scf(3);
clf();
plot(t,resp(4,:),'b-');
plot(t,z(3,:),'r-');
plot(t,X_KF(4,:),'y-');


axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;
xtitle("Kalman - Bank Angle");

//figure 4 - Process Covariance
scf(4);
clf();
plot(t,squeeze(Pr_cov(2,2,:))','y-');
plot(t,squeeze(Pr_cov(2,3,:))','g-');
plot(t,squeeze(Pr_cov(3,3,:))','r-');
plot(t,squeeze(Pr_cov(4,4,:))','b-');

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;
xtitle("Process covariance");
