% Arcticle Title: "A New Iterative Weighted Least Squares Algorithm for 1-D SA Localization"
% Publised on IEEE Signal Processing, 2025
% Codes provided by Yanbin Zou, April, 7, 2025
% Copyright (C) 2025
%       Department of Eletronic Engineering
%       Shantou University 
%       Shantou, Guangdong, 515063, China.
%       ybzou@stu.edu.cn
% The code can reproduce the result in Fig.5 of the reference paper.
NN2=6;
N=1000;%%%Number of Monte Carlo
 S=[ 0   500 -500 500  -500
     0   500  500 -500 -500
     0    0   0   0    0]; %%% Linear arrays' locations
ag=[-2.61,2.78,-0.54,-0.17,-2.97]; %%% Angles of attitude vector 
ae=[0.92,0.19, -0.64, 0.16, -0.73];
M=size(S,2);G=zeros(3,M);
for i=1:M
    G(:,i)=[cos(ae(i))*cos(ag(i));cos(ae(i))*sin(ag(i));sin(ae(i))];%%%Attitude vector 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U0=zeros(3,N);%%%%% Random target locations
for i=1:N
    U0(:,i)=1000*2*(rand(3,1)-0.5*ones(3,1));
end
%%%%%%%%%%%%%%%%%%%%%%%
for K=1:NN2
    K
    k=10^((K-8)/2);    
    Qe=k^2*eye(M);Noise1=zeros(N,M);
    for i=1:N
        Noise1(i,:)=mvnrnd(zeros(1,M),Qe);
    end   
    MSE=0;MSE1=0;         
        for n=1:N
            u0=U0(:,n);
            psi=zeros(M,1);           
            for i=1:M
                psi(i)=acos( (u0-S(:,i))'*G(:,i)/norm(u0-S(:,i)) )+Noise1(n,i);
            end
            %%%%%CRLB
            q=zeros(3,M);
            for i=1:M
                q(:,i)=-1/sqrt(1-((u0-S(:,i))'*G(:,i)/norm(u0-S(:,i)))^2)*...
                    (  norm(u0-S(:,i))*G(:,i)-(u0-S(:,i))*(u0-S(:,i))'*G(:,i)/norm(u0-S(:,i)) )/(norm(u0-S(:,i))^2);
            end
            I=q*inv(Qe)*q';
            CRB=inv(I);
            MSE=MSE+(CRB(1,1)+CRB(2,2)+CRB(3,3)  );
            

            %%%%%%%%% Multiple initial values                    
            M1=20;%%%%% Number of initial values
            x0=[];
            for m=1:M1
                xe0=[-1000;-1000;-1000]+2000*rand(3,1); %% Target space of interest 
                x0=[x0,xe0];
            end
            u= IWLS_Zou_SPL(G,S,psi,Qe,x0);
            MSE1=MSE1+norm(u-u0)^2;              
    end
    RMSE(K)=sqrt(MSE/N);
    RMSE1(K)=sqrt(MSE1/N);
end
figure(1)
K=1:NN2;
semilogy((K-8)*5,RMSE,'r-.',(K-8)*5,RMSE1,'b-s','LineWidth',1)
grid on
xlabel('\sigma (dB rad)' );
ylabel('RMSE(m)');
legend('CRLB','Proposed','location','northwest')
ylim([0,1000]);


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



