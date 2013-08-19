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
pnoise = [t;0.01*wgn(num_samples+1,1,0)']';
onoise = [t;0.01*wgn(num_samples+1,1,0)']';
[b,a] = butter(order,fnorm);
sim('butterworth_mdl.mdl')
plot(t,y,t,y_n,t,u(:,2))

e_noise = 0.01*wgn(num_samples+1,1,0);

%MCMC part

D = filter(b,a,u(:,2));
b = [-0.3,0.2];
a = [1,0.6];
y = filter(b,a,u(:,2));

chi1 = sum((D-y).^2);
sigma = 0.24;
a_save(10000,order+1) = 0;
b_save(10000,order+1) = 0;
chi_save(10000) = 0;
ii = 0;

accepted = 0;
for n = 1:90000
    
    a1 = a + sigma*randn(1,order+1);
    b1 = b + sigma*randn(1,order+1);
    a1(1) = 1;
    y_cand = filter(b1,a1,u(:,2));
    
    chi2 = sum((D-y_cand).^2);
    ratio = exp((-chi2+chi1))/(2*0.01^2);
    
    if rand < ratio
        a = a1;
        b = b1;
        chi1 = chi2;
        accepted = accepted+1;
    end
       
    a_save(n,:) = a1;
    b_save(n,:) = b1;
    chi_save(n) = chi1;
    
    
end

min_loc = find(chi_save==min(min(chi_save)));
a2 = a_save(min_loc(1),:);
b2 = b_save(min_loc(1),:);
y = filter(b2,a2,u(:,2));
[b,a] = butter(order,fnorm);
y_ideal = filter(b,a,u(:,2));
figure
plot(t,y,t,y_ideal)
a_ratio = accepted/n
