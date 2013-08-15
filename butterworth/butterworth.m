clear all
order = 1;
fs = 400;
fc = 20;
fnorm = fc*2/fs;
dt = 1/fs;
t1 = 1.0;
t0 = 0;
t = t0:dt:t1;
num_samples = (t1-t0)/dt;
u = [t;ones(num_samples+1,1)']';
u(1,2) = 0;
pnoise = [t;0.05*wgn(num_samples+1,1,0)']';
onoise = [t;0.05*wgn(num_samples+1,1,0)']';
[b,a] = butter(order,fnorm);
sim('butterworth_mdl.mdl')
plot(t,y,t,y_n,t,u(:,2))