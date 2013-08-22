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
%sim('butterworth_mdl.mdl')
%plot(t,y,t,y_n,t,u(:,2))

e_noise = 0.01*wgn(num_samples+1,1,0);

%MCMC part

D = filter(b,a,u(:,2));
b = [0.1367,0.1367];
a = [1.0,rand];
y = filter(b,a,u(:,2));

chi1 = sum((D-y).^2);
sigma = 1.0;
ii = 0;
flg = 0;
accepted = 0;
for n = 1:10
    
    a1 = a + sigma*randn(1,order+1);
    b1 = b + sigma*randn(1,order+1);
    a1(1) = 1;
    %a1(2) = -0.7265;
    b1(1) = 0.1367;
    b1(2) = 0.1367;
    y_cand = filter(b1,a1,u(:,2));
    a1(2)
    b1;
    chi1
    chi2 = sum((D-y_cand).^2)
    ratio = exp(-chi2+chi1);
    
    if rand < ratio
        a = a1;
        b = b1;
        chi1 = chi2;
        accepted = accepted+1;
    end
    
    if(mod(ii,500)==0 && ii ~=0 && flg==0)
        accepted
        if(accepted/ii < 0.3)
            sigma = sigma/3.2;
            ii = 0;
            accepted = 0;
        elseif(accepted/ii > 0.4)
            sigma = sigma*3.2;
            ii = 0;
            accepted = 0;
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
y = filter(b2,a2,u(:,2));
[b,a] = butter(order,fnorm);
y_ideal = filter(b,a,u(:,2));
figure
%plot(t,y,t,y_ideal)
hist(a_save(:,2),100)
a_ratio = accepted/(n-burnin)
sigma
mean(b_save(:,1))
sqrt(var(b_save(:,1)))
