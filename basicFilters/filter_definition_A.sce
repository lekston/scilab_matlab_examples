    //  Second order linear state space system derived from a transfer function
    //      with natural frequency  1Hz = 6.28 rad/s , damping 0.9 , 
    //      sample rate of 50Hz
    //
    //      x_n = A*x + b*u
    //      y   = c*x + d*u
    //

Order   = 2;    // The order of the filter
Fs      = 50;   // The sampling frequency
Fcutoff = 1;    // The cutoff frequency

// We design a low pass elliptic filter
hz = iir(Order,'lp','ellip',[Fcutoff/Fs/2 0],[0.1 0.1]);

// We compute the frequency response of the filter
[frq,repf]=repfreq(hz,0:0.001:0.5);
[db_repf, phi_repf] = dbphi(repf);

// And plot the bode like representation of the digital filter
scf(2);
clf();
subplot(2,1,1);
plot2d(Fs*frq,db_repf);
xtitle('Obtained Frequency Response (Magnitude)');
subplot(2,1,2);
plot2d(Fs*frq,phi_repf);
xtitle('Obtained Frequency Response (Phase in degree)');

ss_iir = tf2ss(hz);

s=poly(0,'s');
tf = (4*%pi*%pi)/(s^2 + 2*0.9*2*%pi*s + 4*%pi*%pi);
ss = tf2ss(tf);
ss_dscr = dscr(ss, 0.02); //50Hz

t = 0:0.02:20;
u = zeros(1,1001);

for i = 1:1001
u(i) = sin(i*%pi/180);// sin(i*%pi/5)*0.2;
end

//flt50hz = dsimul(ss_dscr,u);
flt50hz = flts(u,ss_dscr);

scf(3);
clf();
plot(t,flt50hz(:),'b-');
plot(t,u(:),'g-');
xtitle("filter 50Hz");


axis = gca();
axis.data_bounds = [0 -1.2; 20 1.2];
