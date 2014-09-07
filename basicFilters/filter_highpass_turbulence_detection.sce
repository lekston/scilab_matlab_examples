    //  Second order linear state space system derived from a transfer function
    //      with natural frequency  1Hz = 6.28 rad/s , damping 0.9 , 
    //      sample rate of 50Hz
    //
    //      x_n = A*x + b*u
    //      y   = c*x + d*u
    //

T_0 = 10; //0.1Hz
T_1 = 2; //0.5Hz
T_2 = 0.2; //5Hz
T_3 = 0.05; //20Hz

k_T0 = 1;
k_T1 = 0;
k_T2 = 0.2;
k_T3 = 0.2;

T_missingSample = 2; // in i, not in t(i)

f_dscr = 50; //Hz

t = 0:0.02:20;
t25hz = 0:0.04:20;
u = zeros(1,1001);

//design highpass
//F_s = 0.3 Hz = 1.88 rad/s       Rs = 40 dB
//F_p = 1.0 Hz = 6.28 rad/s        Rp = 0.2 dB
om = [1.88, 6.28];
deltap = 0.2;
deltas = 0.01;

//generate approx - CRASHES SCILAB CONSOLE!!!
[cells,fact,zzeros,zpoles]=eqiir('hp','ellip',om,deltap,deltas);
//h=fact*poly(zzeros,'z')/poly(zpoles,'z');

// Generate a low pass filter
[_zeros,pols,gain] = zpell(3,60,10,50);

s=poly(0,'s');
tf = (4*%pi*%pi)/(s^2 + 2*0.9*2*%pi*s + 4*%pi*%pi);
ss50hz = syslin('c', tf);
ss50hz_dscr = dscr(ss50hz, 1/f_dscr); //0.02

ss_k1 = syslin('c', s/s);
//lowpass filter in negative feedback loop (giving high pass)
ss_hp = ss_k1 - ss50hz;

ss_hp50hz_dscr = dscr(ss_hp, 1/f_dscr); //0.02

for i = 1:1001
u(i) = sin(t(i)*2*%pi/T_0)*k_T0;
u(i) = u(i) + sin(t(i)*2*%pi/T_1)*k_T1;
u(i) = u(i) + sin(t(i)*2*%pi/T_2)*k_T2;
u(i) = u(i) + sin(t(i)*2*%pi/T_3)*k_T3;
end

flt50hz = dsimul(ss50hz_dscr,u);

flt_hp50hz = dsimul(ss_hp50hz_dscr,u);

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


scf(1);
clf();
plot(t,flt50hz(:),'b-');
plot(t,u(:),'g-');
plot(t,flt_hp50hz(:),'r-');

//plot(t,norm_flt50hz(:),'y-');
xtitle("filter 50Hz");


axis = gca();
axis.data_bounds = [0 -2; 20 2];
