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


clear all;

//Generating Noise(non Gaussian)
T_0 = 10; //0.1Hz
T_1 = 2; //0.5Hz
T_2 = 0.2; //5Hz
T_3 = 0.0125; //80Hz

k_T0 = 0.02;
k_T1 = 0.03;
k_T2 = 0.05;
k_T3 = 0.4;

f_dscr = 50; //Hz

t = 0:(1/f_dscr):20;
v = zeros(1,1001);

for i = 1:size(t,2)
v(i) = sin(t(i)*2*%pi/T_0)*k_T0;
v(i) = v(i) + sin(t(i)*2*%pi/T_1)*k_T1;
v(i) = v(i) + sin(t(i)*2*%pi/T_2)*k_T2;
v(i) = v(i) + sin(t(i)*2*%pi/T_3)*k_T3;
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

//Initial conditions: 1.7deg/s yaw rate; 3deg bank
X0 = [0*(%pi/180);
      1.7*(%pi/180);
      0*(%pi/180);
      3*(%pi/180)];
          
SI = syslin('c', A, B, C, D, X0);

u = zeros(2, size(t,2));
for i = 1:size(t,2)
    u(:,i) = [3*(%pi/180); -6*(%pi/180)]; 
//quasi static at xi = 1deg, dzeta = -4.45
end
resp = csim(u,t,SI);



// ************************** KALMAN FILTER ALGORITHM ********************

//initialization
SI_dscr = dscr(SI,(1/f_dscr));

A_dscr = SI_dscr(2);

X_est = X0;
X_KF = X0;
Inn(1) = 0;
S(1) = 0;

//generating measurement noise
//v = grand(1,size(t,2),'unf',-1,1);

//computing measurement noise covariance
[R, MeasMean] = corr(v,1);

//assuming model noise covariance
Q = R/10;

//assuming initial value of the process noise
Pr_cov = eye(4,4)*R/10;
Pr_cov_est = Pr_cov;

//since measurement is limited to bank angle; the measurement matrix H != Id
H = [0 0 0 1];

//calculating filter output
for i = 1:(size(t,2)-1)
    
    //predict the state
    X_est(:,i+1) = A_dscr*X_KF(:,i);
    
    //predict the state covariance
    Pr_cov_est(:,:,i+1) = A_dscr*Pr_cov(i)*A_dscr' + H*Q*H';



    //noisy measurement
    z(i+1) = resp(4,i+1) + v(i+1)';
    
    //Innovation
    Inn(i+1) = z(i+1) - H*X_est(:,i+1);
    
    //Covariance of the innovation
    S(i+1) = R + H*Pr_cov_est(:,:,i+1)*H';
    
    //update Kalman gain
    K(:,i+1) = Pr_cov_est(:,:,i+1)*H'*inv(S(i+1));


    
    //update state
    X_KF(:,i+1) = X_est(:,i+1) + K(:,i+1)*Inn(i+1);
    
    //update the state covariance
    Pr_cov(:,:,i+1) = Pr_cov_est(:,:,i+1) - K(:,i+1)*S(i+1)*(K(:,i+1)')
    
    
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
plot(t,z(:),'r-');
plot(t,X_KF(4,:),'y-');
plot(t,resp(4,:),'b-');
xtitle("Kalman - Bank Angle");

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;

//figure 4 - Process Covariance
scf(4);
clf();
plot(t,squeeze(Pr_cov(2,2,:))','y-');
plot(t,squeeze(Pr_cov(3,3,:))','g-');
plot(t,squeeze(Pr_cov(4,4,:))','b-');
xtitle("Process Covariance");

axis = gca();
axis.foreground = -2;
axis.background = -1;
axis.font_color = -2;
axis.mark_foreground = -2;
axis.mark_background = -1;

