    //  Second order linear state space system derived from a transfer function
    //      with natural frequency  1Hz = 6.28 rad/s , damping 0.9 , 
    //      sample rate of 50Hz
    //
    //      x_n = A*x + b*u
    //      y   = c*x + d*u
    //
    //   Matlab code:
    //      fcm_tf = tf(2*pi*2*pi,[1, 2*0.9*2*pi,2*pi*2*pi])
    //      fcm_ss = ss(fcm_tf)
    //      fcm_dis = ss(c2d(fcm_tf, 1/50,'matched'))
    //      x_init = (eye(size(fcm_dis.a))-fcm_dis.a)\fcm_dis.b


T_0 = 10; //0.1Hz
T_1 = 2; //0.5Hz
T_2 = 0.2; //5Hz
T_3 = 0.0125; //80Hz

k_T0 = 0.5;
k_T1 = 1;
k_T2 = 0.1;
k_T3 = 3;

tau_T0 = 15; //in i, not in t(i)
T_missingSample = 2; // in i, not in t(i)

f_dscr = 50; //Hz

t = 0:0.02:20;
t25hz = 0:0.04:20;
u = zeros(1,1001);
s=poly(0,'s');
tf = (4*%pi*%pi)/(s^2 + 2*0.9*2*%pi*s + 4*%pi*%pi);
ss50hz = syslin('c', tf);
ss50hz_dscr = dscr(ss50hz, 1/f_dscr); //0.02

ss25hz = syslin('c', tf);
ss25hz_dscr = dscr(ss25hz, 2/f_dscr); //0.04


for i = 1:1001
u(i) = sin(t(i)*2*%pi/T_0)*k_T0;
u(i) = u(i) + sin(t(i)*2*%pi/T_1)*k_T1;
u(i) = u(i) + sin(t(i)*2*%pi/T_2)*k_T2;
u(i) = u(i) + sin(t(i)*2*%pi/T_3)*k_T3;
end

flt50hz = dsimul(ss50hz_dscr,u);

//every second sample substituted by zero hold
for i = 1:1001
    if (modulo(i,T_missingSample) == 0) then
        u_zh(i) = u(i-1);
    else
        u_zh(i) = u(i);
    end
end
flt50hz_zh = dsimul(ss50hz_dscr,u_zh');
//flt50hz = flts(u,ss50hz_dscr); //same as above

//only every second sample fed to 25hz filter
for i = 1:501
    u25hz(i) = u(2*i - 1);
end
flt25hz = dsimul(ss25hz_dscr,u25hz');

//normalize by removing lowest frequency:
for i = 1:(1001 - tau_T0)
norm_flt50hz(i+tau_T0) = flt50hz(i+tau_T0) - sin(t(i)*2*%pi/T_0)*k_T0;
end

scf(1);
clf();
plot(t,flt50hz(:),'b-');
plot(t,u(:),'g-');
plot(t,flt50hz_zh(:),'r-');
plot(t25hz,flt25hz(:),'y-');
//plot(t,norm_flt50hz(:),'y-');
xtitle("filter 50Hz");


axis = gca();
axis.data_bounds = [0 -2; 20 2];