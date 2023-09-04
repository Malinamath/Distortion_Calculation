iter=100; zi=0.02;
z=0.01:0.005:0.1;
d= [0.3500    0.3300    0.3400    0.3450    0.3350    0.3250    0.3050    0.3400    0.3350    0.3350    0.3200    0.3700 0.3250    0.3000    0.3000    0.3100    0.3050    0.3650    0.2850];


plot(z,d);
grid on;

%title('2-D Line Plot')
xlabel('zeta')
ylabel('Distortion')
%legend('Shanon Limit','400,100','T2{(70,3,3;g28)}','Time Sharing');
