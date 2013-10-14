N=50;
elim = 0.1;
M=1e5;


u=ones(1,N);
a=[1,-0.9];
b=sum(a);
z = filter(b,a,u);

e = elim*(2*rand(size(z))-1);
y=z+e;

theta = [b,a(2)]';

THETA = zeros(2,M);  THETA(:,1)=theta;

% Generate a proposal;

for k=2:M

  xi = THETA(:,k-1) + 0.001*randn(size(theta));
  ytest = filter(xi(1),[1,xi(2)],u);

  if max(abs(y-ytest))>elim
      THETA(:,k) = THETA(:,k-1);
  else
      THETA(:,k) = xi;
  end;


end;

plot(THETA(2,:))

figure(2)
hist(THETA(1,:),100);
figure(3)
hist(THETA(2,:),100);

