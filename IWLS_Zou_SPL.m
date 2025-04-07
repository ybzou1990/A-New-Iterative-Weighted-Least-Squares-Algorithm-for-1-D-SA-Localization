function [u] = IWLS_Zou_SPL(G,S,psi,Qe,x0)
%%  G: A fat matrix constituted by attitude vectors
%%% S: Linear arrays' location
%%% psi: 1-D space angle 
%%% Qe: Noise covariance matrix
%%% x0: A fat matrix constituted by multiple iniital values
%%%%%%%   Propsosed  IWLS
          M=size(S,2);%%% Number of linear arrays
          M1=size(x0,2);%%% Number of initial values
          cf=[];X=[];%%%           
            for m0=1:M1
                u=x0(:,m0);
                A=zeros(M,3);
                b=zeros(M,1); threshold=0.01; W=eye(M);
                for m=1:10
                    ua=u;
                    for i=1:M
                        A(i,:)=(cos(psi(i))*(u-S(:,i))/norm(u-S(:,i))-G(:,i))';
                        b(i)=(cos(psi(i))*(u-S(:,i))/norm(u-S(:,i))-G(:,i))'*S(:,i);
                    end                    
                    u=inv(A'*W*A)*A'*W*b;
                    if norm(u-ua)<threshold
                        break
                    else                        
                        D=zeros(M,M);
                        for i=1:M
                            D(i,i)=norm(u-S(:,i))*sin(psi(i));
                        end
                        W=inv(D*Qe*D');                        
                    end                    
                end
                d0=zeros(M,1);
                for i=1:M
                    d0(i)=G(:,i)'*(u-S(:,i))/norm(S(:,i)-u)-cos(psi(i));
                end
                cf(m0)=d0'*d0;
                X=[X,u];                
            end
            [f1,ind]=min(cf);
            u=X(:,ind); %%% estimated target location          
end

