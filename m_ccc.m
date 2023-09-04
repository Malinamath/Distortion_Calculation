function b=m_ccc(m,a,u)
[p q]=size(a);x=max(max(a));num=1;
b=sparse(x*m,p*m);
for i=1:p
    
        if a(i,1)~=0            
            for k=1:m
                
                    b((a(i,1)-1)*m+k,(i-1)*m+k)=1;
                
            end
        end
    
    for j=2:length(a(i,:))
        if a(i,j)~=0
            s=u(num);num=num+1;
            for k=1:m
                if k+s<=m
                    b((a(i,j)-1)*m+k+s,(i-1)*m+k)=1;
                else
                     b((a(i,j)-2)*m+k+s,(i-1)*m+k)=1;
                end
            end
        end
    end
end