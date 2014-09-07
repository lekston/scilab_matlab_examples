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

s=poly(0,'s');
tf = (4*%pi*%pi)/(s^2 + 2*0.9*2*%pi*s + 4*%pi*%pi);
ss = syslin('c', tf);
ss_dscr = dscr(ss, 0.02); //50Hz

t = 0:0.02:20;
u = zeros(1,1001);

T_0 = 10; //0.1Hz
T_1 = 2; //0.5Hz
T_2 = 0.2; //5Hz
T_3 = 0.05; //20Hz

k_T0 = 0.5;
k_T1 = 1;
k_T2 = 0.1;
k_T3 = 0.1;

tau_T0 = 15; //in i, not in t(i)
T_missingSample = 6; // in i, not in t(i)

for i = 1:1001
u(i) = sin(t(i)*2*%pi/T_0)*k_T0;
u(i) = u(i) + sin(t(i)*2*%pi/T_1)*k_T1;
u(i) = u(i) + sin(t(i)*2*%pi/T_2)*k_T2;
u(i) = u(i) + sin(t(i)*2*%pi/T_3)*k_T3;
end

flt50hz = dsimul(ss_dscr,u);

//every second sample substituted by zero hold
for i = 1:1001
    if (modulo(i,T_missingSample) == 0) then
        u_zh(i) = u(i-1);
    else
        u_zh(i) = u(i);
    end
end

flt50hz_zh = dsimul(ss_dscr,u_zh');

//flt50hz = flts(u,ss_dscr); //same as above

for i = 1:986
norm_flt50hz(i+tau_T0) = flt50hz(i+tau_T0) - sin(t(i)*2*%pi/T_0)*k_T0;
end

scf(1);
clf();
plot(t,flt50hz(:),'b-');
//plot(t,u(:),'g-');
plot(t,flt50hz_zh(:),'r-');
plot(t,norm_flt50hz(:),'y-');
xtitle("filter 50Hz");


axis = gca();
axis.data_bounds = [0 -1.2; 20 1.2];