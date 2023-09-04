function [total_average] = traialDistortion(zi,iter,G,S,trial)
%H = alist2sparse('za21') ;
%G=H';
arr_distortion=[];
average=0;
for tri=1:trial
    distortion=Distortion(zi,iter,G,S);
    tri=tri+1;
average=average+distortion;
arr_distortion=[arr_distortion,distortion];
end 
total_average=average/trial;
disp(arr_distortion);
end