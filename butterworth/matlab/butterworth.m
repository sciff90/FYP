clear all
order = 1;
fs = 400;
fc = 20;
fnorm = fc*2/fs;
dt = 1/fs;
t1 = 0.15;
t0 = 0;
t = t0:dt:t1;
num_samples = (t1-t0)/dt;
u = ones(num_samples+1,1)';
u(1) = 0;
u = [ones(num_samples+1,1)';3.*ones(num_samples+1,1)';-2.*ones(num_samples+1,1)';cos(-10.*t)];
u(:,1) = 0;
u = u';
[b,a] = butter(order,fnorm);

%MCMC part

D = filter(b,a,u);
b_curr = rand(1,order+1);
a_curr = rand(1,order+1);
a_curr(1) = 1;
a_best = a_curr;
b_best = b_curr;
y_curr = filter(b_curr,a_curr,u);

chi1 = norm(sum((D-y_curr).^2));
chi_best = chi1;
sigma = 1.0;
ii = 0;
flg = 0;
accepted = 0;
n = 0;

while(n<=80000)
    
    a_cand = a_curr + sigma*randn(1,order+1);
    b_cand = b_curr + sigma*randn(1,order+1);
    a_cand(1) = 1;
    %a_cand(2) = a(2);
    %a_cand(3) = a(3);
    %b_cand(1) = b(1);
    %b_cand(2) = b(2);
    %b_cand(3) = b(3);
    y_cand = filter(b_cand,a_cand,u);
    chi1;
    chi2 = norm(sum((D-y_cand).^2));
    ratio = exp(-(chi2)+(chi1));
    
    if rand < ratio
        a_curr = a_cand;
        b_curr = b_cand;
        chi1 = chi2;
        if(chi2<chi_best)
            chi_best = chi2;
            a_best = a_curr;
            b_best = b_curr;            
        end
        accepted = accepted+1;
    end
    
    if(mod(ii,1000)==0 && ii ~=0 && flg==0)
        accepted;
        if(accepted/ii < 0.3)
            sigma = sigma/1.2;            
            ii = 0;
            accepted = 0;
        elseif(accepted/ii > 0.4)
            sigma = sigma*1.2;            
            ii = 0;
            accepted = 0;
        else
            burnin = n-1;
            flg = 1;
        end
    end
    ii = ii+1;
    n = n+1;
    if(flg==1)
        a_save(n-burnin,:) = a_curr;
        b_save(n-burnin,:) = b_curr;
        chi_save(n-burnin) = chi1;
    end   

end
format long
mean_err = mean(chi_save)
best_err = chi_best
a_ratio = accepted/(n-burnin)
burnin
sigma
a;
a_best;
a_mean = mean(a_save);
b;
b_best;
b_mean = mean(b_save);

y = filter(b_best,a_best,u);
y_ideal = filter(b,a,u);
figure(1);

subplot(3,2,1);
plot(t,y,t,y_ideal)
title('Output Response')
xlabel('time (s)')

subplot(3,2,2);
hist(a_save(:,2),100)
title(strcat('a[2] mean = ',num2str(mean(a_save(:,2)))))
ylabel('frequency')

subplot(3,2,3);
hist(b_save(:,1),100)
title(strcat('b[1] mean = ',num2str(mean(b_save(:,1)))))
ylabel('frequency')

subplot(3,2,4);
hist(b_save(:,2),100)
title(strcat('b[2] mean = ',num2str(mean(b_save(:,2)))))
ylabel('frequency')

h = axes('Position',[0 -0.4 1 1],'Visible','off');
str(1) = {'Expected Values for a and b Vs best found'};
str(2) = {''};
str(3) = {strcat('Expected a values     = ',num2str(a))};
str(4) = {strcat('Best Found a values = ',num2str(a_best))};
str(5) = {strcat('Expected b values     = ',num2str(b))};
str(6) = {strcat('Best Found b values = ',num2str(b_best))};
set(gcf,'CurrentAxes',h)
text(.025,.6,str,'FontSize',18)
 