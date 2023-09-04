

function f = ziDistortion(x)
iter=100;
m=5;
a=[1 4 0 0 0 0 0 0 ;2 9 0 0 0 0 0 0 ;3 10 0 0 0 0 0 0 ;5 6 0 0 0 0 0 0 ;7 8 0 0 0 0 0 0 ;1 5 0 0 0 0 0 0 ;2 6 8 0 0 0 0 0 ;3 4 9 0 0 0 0 0 ;1 7 10 0 0 0, 0 0 ;2 5 10 0 0 0 0 0 ;3 6 7 0 0 0 0 0 ;4 8 9 10 0 0 0 0 ;1 3 8 9 0 0 0 0 ;2 3 4 5 6 7 9 10 ;1 2 4 5 6 7 8 10 ;1 2 3 4 6 7 8 9 ;2 3 4 5 6 7 9 10 ;1 2 3 4 5 8 9 10 ;1 4 5 6 7 8 9 10 ;1 2 3 5 6 7 8 0];
u=[0,3,4,4,1,0,3,1,0,3,2,0,3,1,3,2,0,2,3,3,4,3,0,2,0,2,2,0,1,0,1,2,4,2,4,1,2,2,1,0,3,2,4,2,2,3,4,2,0,1,2,2,4,2,2,3,0,1,4,0,4,1,2,2,3,1,3,0,1,0];
h=m_ccc(m,a,u);
G=h';
p=0.5;
row=100;
%S=binornd(1,p*ones(row,1)); % generate Bernoulli Process p=0.5
S=[0     0     1     0     1     1     1     1     0     0     0     0     1     0     0     0     0     1     1     0     1     1     1   0     1     1     0     0     1     0     1     1     0     0     0     1     0     1     0     0     1     0     0     0     1     0  1     1     1     0     1     0     0     0     1     1     1     0     1     1     1     0     0     1     1     0     1     1     0   0     0     0     1     0     0     0     0     1     0     1     1     0     1     1     1     1     0     0     1     1     1     0   0     0     1     1     1     0     1     1];

%disp(S);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zhat=ziDecimation(x,G,S,iter);

function [zhat]=ziDecimation(x,G,S,iter)

INDX=0;
dim1=size(G);
rows1=dim1(1);
cols1=dim1(2);
zhat(1,1:cols1)=0;
zero(1,1:rows1)=0;
%V=1:cols1;

s=struct('Bmn',1,'qmn',0,'qm',0,'rmn',1,'Bm',1);

%associate this structure with all non zero elements of G

%initialization
for i= 1:rows1
   for j=1:cols1
      newh(i,j)=s;
      if G(i,j)==1
         newh(i,j).Bmn=0.0500;
         newh(i,j).qmn=0.1;
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%STARTING ITERATION%%%%%%%%%%%%%%%%%%%%%%%%%%55555
for iteration=1:iter
 
    
if(INDX~= length(zhat))
   
    dim=size(G);
rows=dim(1);
cols=dim(2);
%gama(1,1:cols)=1;
%beta=tanh(gama);
mu=1./x;
beta=(1-x)./(1+x);

%initialization for iteration 2
if (iteration==2)
for i= 1:rows
   for j=1:cols
     
      if G(i,j)==1
         newh(i,j).qmn=0;
        % disp(newh(i,j).qmn);
      end
   end
end
end  

%%%%%%%%%%%%%%%%%%%%horizontal step%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i1=1:rows
   	ones_in_row=0;
   	index=1;
   	for j1=1:cols  %finding the indices of ones in each rows(return column indeces)
      	if G(i1,j1)==1
         	ones_in_row(index)=j1;
         	index=index+1;
         	newh(i1,j1).fmn=(1./mu).* newh(i1,j1).qmn;
         %  fprintf('newh(i1,j1).fmn = %d\n',newh(i1,j1).fmn);
           
        end % for if
        
    end %for j1
    
   	for j1=1:index-1
        	tmn=1;
      	for k=1:index-1
         	if k~=j1  %exclude index
            	tmn= tmn .* newh(i1,ones_in_row(k)).Bmn;
                
            end %for if
            
        end %for k
         TMN=(2.*(-1).^((S(i1)+index-1))).* atanh(beta.*tmn);
          
newh(i1,ones_in_row(j1)).rmn=newh(i1,ones_in_row(j1)).fmn+TMN;
  % fprintf('newh(i,j).rmn = %d\n',newh(i1,ones_in_row(j1)).rmn);
    end %for j1
 
end %for i1
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%vertical step%%%%%%%%%%%%
      cols_with_all_zeros = find(all(G==0));
    maxvalue=0;
     Z=[];
	for j=1:cols 
       if(~ismember(j,cols_with_all_zeros))  
   	ones_in_col=0;
   	index=1;
   	for i=1:rows
      	if G(i,j)==1
         	ones_in_col(index)=i; %finding the indices of one in columns(return row indices)
         	index=index+1;
        
      	end
    end
   
   	for i=1:index-1
      	SUM_rmn=0;
        SUM_rm=0;
             
       
      	for k=1:index-1
                SUM_rm=SUM_rm +newh(ones_in_col(k),j).rmn;
         	if k~=i
            	SUM_rmn=SUM_rmn+newh(ones_in_col(k),j).rmn; 
              
            end 
            
        end 
        newh(ones_in_col(i),j).qm=SUM_rm;
        newh(ones_in_col(i),j).qmn=SUM_rmn;
        newh(ones_in_col(i),j).Bmn=tanh(0.5.*SUM_rmn);
        newh(ones_in_col(i),j).Bm=tanh(0.5.*SUM_rm);
        
       
        
  %    fprintf('newh(i,j).qm = %d\n',newh(ones_in_col(i),j).qm);
 %  fprintf('newh(i,j).qmn = %d\n',newh(ones_in_col(i),j).qmn); 
   %  fprintf('newh(i,j).Bmn = %d\n',newh(ones_in_col(i),j).Bmn);
   % fprintf('newh(i,j).Bm = %d\n',newh(ones_in_col(i),j).Bm);
    
 if(abs(newh(ones_in_col(i),j).Bm )> maxvalue) 
         maxvalue=abs(newh(ones_in_col(i),j).Bm); 
         maxx= newh(ones_in_col(i),j).Bm;
       y = j;
       Z = [Z,y];
      
  %disp(Z); 
  

 end
   
    end  %for i
    
   end %if not column zero  
   moxim=maxx;
    end %for j
   
  if(length(Z)>1)
     I= randsample(Z,1);
    % fprintf('I is=%d\n',I); 
  else
      I=Z;
     % fprintf('I is=%d\n',I); 
  end    
 
[INDX,zhat] = decimatezhat(I,zhat,moxim,cols); 
%disp(V);
% disp(zhat);
  S=sourceupdate(S,G,I,zhat);
%disp(S);
 % fprintf('S is=%d\n',S);
   G = graphdecimation(G,I);
 %disp(G);
 
    iteration;
    
end  %end if
end % end for

end  % end main function



%+++++++++++++++++++++++++++++
 function [INDX,zhat] = decimatezhat(I,zhat,moxim,cols)
%cols=4;
%I=[2,4];
%maxvalue=0.3;
%zhat(1,1:cols)=0;

INDX=0;
 for j=1:cols
         if (ismember(j, I))
        
             if  moxim>0 
       zhat(j)=1;
      % fprintf('hello world \n');
       INDX=INDX+1;
             elseif  moxim<0
         zhat(j)=0;
        % fprintf('hello mom \n');
         INDX=INDX+1;
            end %END IF
            
        
         end %first if
 end  %for
 end% for function
%+++++++++++++++++++++++++++++++++++++++++++

  function [S]=sourceupdate(S,G,I,zhat)

%G=[1 0 1 1; 0 1 1 0; 1 0 1 0; 1 1 1 0];
%S=[1 1 0 1];
%I=[2];
%zhat=[0 1 0 1 ];
dim=size(G);
rows=dim(1);
cols=dim(2);

Shat(1,1:rows)=0;
%zero(1,1:rows)=0;

for j=1:cols
   
 if (ismember(j, I))
     for i=1:rows
         
     if G(i,j)==1
         index=1;
         for k=1:cols
             
             if G(i,k)==1
            ones_in_row(index)=k;
         	index=index+1;
           
             end
             
         end
        %  fprintf('index is %d\n',index);
        % if( index>=3)
             Shat(i)=xor(S(i),zhat(j)) ;
             S(i)=Shat(i);
            % fprintf(' Shat(i) is %d\n', Shat(i));
            % fprintf(' S is %d\n', S);
             
         % changem(S,Shat(i),S(i));
          
       %  end
    end  
     end
     
 end
 
end

 end
  %+++++++++++++++++++++++++++
   
function [G] = graphdecimation(G,I)
G(:,I)=0; 
  % disp(G);
  end
   %++++++++++++++++++++++++++
%%%********************************************************************8

%disp(zhat);
zt=(zhat)';
SBAR = mod(G*zt,2);
%fprintf('Shat is=****************\n');
%disp(SBAR);

sum=0;
for i=1:length(S)
    
        if S(i)==SBAR(i)
            diff=0;
        else 
             diff=1;
        end

        sum=sum+diff;
   
end

f=0.5*(1/length(S))*sum;

fprintf('AVERAGE DISTORTION IS=%d\n',f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
