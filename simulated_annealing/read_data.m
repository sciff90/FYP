fid = fopen('Data/climate.dat');
A = fscanf(fid, '%g %g %g', [3 inf]);
fclose(fid);
A = A'
t = A(:,1);
u_in = A(:,2);
y_out = A(:,3);

plot(t,y_out,t,u_in)