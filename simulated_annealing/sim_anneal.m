clear all
order = 2;
fs = 400;
fc = 20;
fnorm = fc*2/fs;
dt = 1/fs;
t1 = 0.15;
t0 = 0;
t = t0:dt:t1;
num_samples = (t1-t0)/dt;
u = [ones(num_samples+1,1)';sin(10.*t)]
stop
u(1,:) = 0;
pnoise = [t;0.01*wgn(num_samples+1,1,0)']';
onoise = [t;0.01*wgn(num_samples+1,1,0)']';
[b,a] = butter(order,fnorm);
%sim('butterworth_mdl.mdl')
%plot(t,y,t,y_n,t,u(:,2))

e1_noise = 0.1*wgn(num_samples+1,1,0);
e2_noise = 0.1*wgn(num_samples+1,1,0);

%MCMC part

D = filter(b,a,u);
b = rand(1,order+1);
a = rand(1,order+1);
a(1) = 1;
y = filter(b,a,u);

chi1 = sum((D-y).^2);
%chi1 = sum(abs(D-y));
sigma = 1.0;
ii = 0;
flg = 0;
accepted = 0;
T1 = 100;
T = T1;
for n = 1:80000
    
    a1 = a + sigma*randn(1,order+1);
    b1 = b + sigma*randn(1,order+1);
    a1(1) = 1;
    %a1(2) = -0.7265;
    %b1(1) = 0.1367;
    %b1(2) = 0.1367;
    y_cand = filter(b1,a1,u);
    a1(2);
    b1;
    chi1;
    chi2 = sum((D-y_cand).^2);
    %chi2 = sum(abs(D-y_cand));
    ratio = exp(-chi2+chi1);
    
    if(rand < exp((chi1-chi2)/T))
        a = a1;
        b = b1;
        chi1 = chi2;
        accepted = accepted+1;
        T = T*0.5;
    end
    
    if(mod(ii,100)==0 && ii ~=0 && flg==0)
        accepted;
        if(accepted/ii < 0.3)
            sigma = sigma/1.2;
            ii = 0;
            accepted = 0;
            T = T1;
        elseif(accepted/ii > 0.4)
            sigma = sigma*1.2;
            ii = 0;
            accepted = 0;
            T = T1;
        else
            burnin = n-1;
            flg = 1;
        end
    end
    ii = ii+1;
    if(flg==1)
        a_save(n-burnin,:) = a1;
        b_save(n-burnin,:) = b1;
        chi_save(n-burnin) = chi1;
    end
    
    
end

accepted
min_loc = find(chi_save==min(min(chi_save)));
a2 = a_save(min_loc(1),:);
b2 = b_save(min_loc(1),:);
y = filter(b2,a2,u);
[b,a] = butter(order,fnorm);
y_ideal = filter(b,a,u);
figure(1);
plot(t,y,t,y_ideal)
%hist(a_save(:,2),100)
%figure
%hist(b_save(:,1),100)
%figure
%hist(b_save(:,2),100)
a_ratio = accepted/(n-burnin)
sigma
mu = mean(b_save(:,1))
sig = sqrt(var(b_save(:,1)))
T
a2 = mean(a_save);
b2 = mean(b_save);
y = filter(b2,a2,u);
[b,a] = butter(order,fnorm);
y_ideal = filter(b,a,u);
figure(2);
plot(t,y,t,y_ideal)
