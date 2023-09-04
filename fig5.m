iter=400; zi=0.02;
x=[100  200 1000  2000 5000 ];
y=[0.2850  0.2825 0.2735 0.2533 0.2524 ];


plot(x,y);
grid on;

%title('2-D Line Plot')
xlabel('Codeword Length')
ylabel('Distortion')
%legend('Shanon Limit','400,100','T2{(70,3,3;g28)}','Time Sharing');
