function [DH] = HammingDistortion(S,Shat)
%The binary Hamming distortion function
for i=1:length(S)
    for j=1:length(Shat)
        if S(i)==Shat(j)
            DH(i,j)=0;
        else 
             DH(i,j)=1;
        end
    end
end

end

