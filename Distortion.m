function [Dbar] = Distortion(zi,iter,G,S)
%p=0.2;
%S=binornd(1,p*ones(row,1)); % generate Bernoulli Process p=0.2

zhat=Decimation(zi,iter,G,S);
disp(zhat);
zt=(zhat)';
SBAR = mod(G*zt,2);
fprintf('Shat is=****************\n');
disp(SBAR);


sum=0;
for i=1:length(S)
    
        if S(i)==SBAR(i)
            diff=0;
        else 
             diff=1;
        end

        sum=sum+diff;
   
end

Dbar=0.5*(1/length(S))*sum;

fprintf('AVERAGE DISTORTION IS=%d\n',Dbar);
end

