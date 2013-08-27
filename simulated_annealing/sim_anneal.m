clear all
order = 2;
fs = 400;
fc = 20;
fnorm = fc*2/fs;
dt = 1/fs;
t1 = 0.1;
t0 = 0;
t = t0:dt:t1;
num_samples = (t1-t0)/dt;
u = [ones(num_samples+1,1)';3.*ones(num_samples+1,1)';-2.*ones(num_samples+1,1)';cos(-10.*t)];
%u = [ones(num_samples+1,1)'];
%u = [ones(num_samples+1,1)'; -1*ones(num_samples+1,1)'];
u(:,1) = 0;
u = u';
pnoise = [t;0.01*wgn(num_samples+1,1,0)']';
onoise = [t;0.01*wgn(num_samples+1,1,0)']';
[b,a] = butter(order,fnorm);
%sim('butterworth_mdl.mdl')
%plot(t,y,t,y_n,t,u(:,2))

d = 0.1;
noise = d*wgn(num_samples+1,4,0);

%MCMC part

D = filter(b,a,u);%+noise;
b_curr = rand(1,order+1);
a_curr = rand(1,order+1);
a_curr(1) = a(1);
d_curr = rand;
d_curr = 0.1;
d_best = d_curr;
y_curr = filter(b_curr,a_curr,u);%+d_curr*wgn(num_samples+1,4,0);

chi1 = sum((D-y_curr).^2);
sigma = 0.001;
ii = 0;
N_iter = 0;
flg = 0;
accepted = 0;
T_max = 100000;
T_change = 0.9999;
T = T_max;
chi_best = chi1;
count = 0;
while(flg==0)
    
    a_cand = a_curr + sigma*randn(1,order+1);
    b_cand = b_curr + sigma*randn(1,order+1);
    temp = d_curr+sigma*(randn(1,1));
    while(temp<0)
        temp = d_curr+sigma*(randn(1,1));
    end
    d_cand = temp;
    d_cand = 0.1;
            
    a_cand(1) = 1;    
    
    y_cand = filter(b_cand,a_cand,u);%+d_cand*wgn(num_samples+1,4,0);    
    
    chi2 = sum((D-y_cand).^2);
        
    if(norm(chi2)<= norm(chi1))
        a_curr = a_cand;
        b_curr = b_cand; 
        d_curr = d_cand;
        if(norm(chi2)<=norm(chi_best))
            a_best = a_curr;
            b_best = b_curr;
            d_best = d_curr;
            chi_best = chi2;
           
        end
        chi1 = chi2;
        accepted = accepted+1;
        T = T*T_change;
    elseif(exp((chi1-chi2)/T)>rand())
        a_curr = a_cand;
        b_curr = b_cand;
        d_curr = d_cand;
        chi1 = chi2;
        accepted = accepted+1;
        T = T*T_change;
    end
       
    %if(mod(ii,10)==0 && ii ~=0 && flg==0)
    %        sigma = sigma/2.5;
    %        if(sigma < 10^-4);
    %            sigma = 10;                
    %        end
    %    ii=0;
    %    accepted = 0;        
    %end
        
    if(norm(chi_best)<10^-3)
            burnin = N_iter;
            flg = 1;
    end
    ii = ii+1;
    N_iter = N_iter+1;
    
    
    if(T<1*10^-9)
        T = T_max;
        b_curr = rand(1,order+1);
        a_curr = rand(1,order+1);
        a_curr(1) = a(1);
        y_curr = filter(b_curr,a_curr,u);
        chi1 = sum((D-y_curr).^2);
        sigma = 0.001;
    end
    norm(chi_best)
end

y = filter(b_best,a_best,u);% + d_best*wgn(num_samples+1,4,0);
[b,a] = butter(order,fnorm);
y_ideal = filter(b,a,u);%+noise;
figure(1);
plot(t,y,t,y_ideal)
format long
a
a_best
b
b_best
%d
%d_best
norm_err = norm(chi_best)
