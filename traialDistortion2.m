%function [total_average] = traialDistortion2(zi,iter,trial,P1)
 function [total_average] = traialDistortion2(zi,iter,trial,N1, M1, Dv, DoYouWantACE)
G = ProgressiveEdgeGrowthV2(N1, M1, Dv, DoYouWantACE);
%H=P1;
%G=H';
arr_distortion=[];
average=0;
for tri=1:trial
    p=0.5;
S=binornd(1,p*ones(100,1)); % generate Bernoulli Process p=0.2
    distortion=Distortion(zi,iter,G,S);
    tri=tri+1;
average=average+distortion;
arr_distortion=[arr_distortion,distortion];
end 
total_average=average/trial;
%disp(arr_distortion);
end