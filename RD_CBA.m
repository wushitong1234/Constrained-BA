function [time_avg_ABP,time_avg_BA,time_avg_AS,Err_1,Err_2,Err_3,num_1,num_2,num_3,num_4,R,zeta] = RD_CBA(D,number,tol,opt)

M=100;
N=M;K=M;snr=1;h=8;
delta=2*h/N;
x=((-h+delta/2):delta:(h-delta/2))';
y=((-h+delta/2):delta:(h-delta/2))';
r0=1/N*ones(N,1);

if (opt.parp==1)
    lambda = opt.lambda;
    p=delta*1/(2*lambda)*exp(-abs(x)/lambda);
    p=p/sum(p);
    d=abs(x*ones(1,N)-ones(M,1)*y');
else
    sigma = opt.sigma;
    p=delta*1/sqrt(2*pi*sigma^2)*exp(-x.^2/(2*sigma^2));
    p=p/sum(p);
    d=(x*ones(1,N)-ones(M,1)*y').^2;
end


d2=d.*d; 


time_tol_ABP=0;
time_tol_BA=0;
time_tol_AS=0;



for i1=1:number  
%%======= CBA
time_ABP=tic;
zeta=1;
w1=1/N*ones(M,N);r1=r0;
R=-1;err=1;num_abp=0;
expzd=exp(-zeta*d);
while (err>tol)
    num_abp=num_abp+1;
    i_abp=0;
    %G=sum(p.*(((exp(-zeta*d).*d)*r1)./(exp(-zeta*d)*r1)))-D;
    %expzd=exp(-zeta*d);
    G=p'*(((expzd.*d)*r1)./(expzd*r1))-D;
    while ((i_abp<30)&&(abs(G)>1e-12))
        i_abp=i_abp+1;
        %expzd=exp(-zeta*d);
        expzdr=expzd*r1;
        %zeta=zeta-(G)/(sum(p.*((((exp(-zeta*d).*d)*r1).^2-((exp(-zeta*d)*r1).*((exp(-zeta*d).*d.*d)*r1)))./((exp(-zeta*d)*r1).^2))));
        zeta=zeta-(G)/(p'*((((expzd.*d)*r1).^2-(expzdr.*((expzd.*d2)*r1)))./(expzdr.^2)));
        %G=sum(p.*(((exp(-zeta*d).*d)*r1)./(exp(-zeta*d)*r1)))-D;
        expzd=exp(-zeta*d);
        G=p'*(((expzd.*d)*r1)./(expzd*r1))-D;
    end
    %w1=(r1').*(exp(-zeta*d));
    w1=(r1').*expzd;
    w1=w1./sum(w1,2);
    %w1=exp(-zeta*d)*diag(r1);
    %w11=diag(1./sum(w1,2));
    %w1=w11*w1;
    r1=w1'*p;
    w1(w1<1e-60)=1e-60; 
    r1(r1<1e-60)=1e-60;
    %R1=sum(sum((diag(p)*w1).*log(w1)))-sum(r1.*log(r1));
    %R1=sum(sum((w1.*p).*log(w1)))-sum(r1.*log(r1));
    R1=sum(p'*(w1.*log(w1)))-sum(r1.*log(r1));
    err=abs(R-R1);R=R1;
end
time_ABP=toc(time_ABP);
time_tol_ABP=time_tol_ABP+time_ABP;
result_ABP=R; 


%%======= BA
time_BA=tic;
b1=0;b2=50;num_ba_out=0;tols=tol/1000;
for j_ba=1:1000
    if abs(b1-b2)<tols
        break;
    end
    zeta=(b1+b2)/2;
    num_ba_out=num_ba_out+1;
    w1=1/N*ones(M,N);r1=r0;
    R=-1;err=1;
    num_ba_in=0;
    while (err>tol)
        num_ba_in=num_ba_in+1;
        w1=(r1').*(exp(-zeta*d));
        w1=w1./sum(w1,2);
        %w1=exp(-zeta*d)*diag(r1);
        %w11=diag(1./sum(w1,2));
        %w1=w11*w1;
        r1=w1'*p;
        w1(w1<1e-60)=1e-60;r1(r1<1e-60)=1e-60;
        %R1=sum(sum((diag(p)*w1).*log(w1)))-sum(r1.*log(r1));
        %R1=sum(sum((w1.*p).*log(w1)))-sum(r1.*log(r1));
        R1=sum(p'*(w1.*log(w1)))-sum(r1.*log(r1));
        err=abs(R-R1);R=R1;
    end
    if(sum(sum(diag(p)*(w1.*d)))<D)
        b2=zeta;
    else
        b1=zeta;
    end
end
time_BA=toc(time_BA);
time_tol_BA=time_tol_BA+time_BA;
result_BA=R; 


%%======= AS
time_AS=tic;
zeta=1;
phi=ones(M,1);r=r0;
K=exp(-zeta*d);R=-1;num_as=0;err=1;
% while (err>tol || num_as<num_abp)
while (err>tol)
    num_as=num_as+1;
    psi=1./(K'*(phi.*p));
    phi=1./(K*(psi.*r));
    i_as=0;
    %G=sum(sum((exp(-zeta*d).*d).*((p.*phi)*(psi.*r)')))-D;
    G=(p.*phi)'*(K.*d)*(psi.*r)-D;
    while ((i_as<30)&&(abs(G)>1e-12))
        i_as=i_as+1;
        %zeta=zeta+(G)/(sum(sum((exp(-zeta*d).*d.*d).*((p.*phi)*(psi.*r)'))));
        zeta=zeta+(G)/((p.*phi)'*(K.*d2)*(psi.*r));
        %G=sum(sum((exp(-zeta*d).*d).*((p.*phi)*(psi.*r)')))-D;
        K=exp(-zeta*d);
        G=(p.*phi)'*(K.*d)*(psi.*r)-D;
    end
    %K=exp(-zeta*d);
    %w=diag(phi)*K*diag(psi.*r);
    %w=((K.*phi)'.*(psi.*r))';
    w=(psi.*r)'.*(K.*phi);
    i_as=0;
    rr=w'*p;
    beta=-log(psi)-0.5;
    eta=max(beta)+1e-5;
    F=sum(rr./(eta-beta))-1;
    while ((i_as<40)&&(abs(F)>1e-14))
        i_as=i_as+1;
        eta=eta+(F)/(sum(rr./((eta-beta).^2)));
        F=sum(rr./(eta-beta))-1;
    end
    r=rr./(eta-beta);
    %ww=diag(phi)*K*diag(psi);
    %ww=((K.*phi)'.*psi)';
    ww=(psi)'.*(K.*phi);
    ww(ww<1e-6)=1e-60;r(r<1e-60)=1e-60;
    % R1=sum(sum((diag(p.*phi)*K*diag(psi.*r)).*log(ww)));
    R1=(p.*phi)'*(K.*log(ww))*(psi.*r);
    err=abs(R-R1);R=R1;
end
time_AS=toc(time_AS);
time_tol_AS=time_tol_AS+time_AS;
result_AS=R; 
end 

time_avg_ABP=time_tol_ABP/number;
time_avg_BA=time_tol_BA/number;
time_avg_AS=time_tol_AS/number;
Err_1=abs(result_ABP-result_BA);
Err_2=abs(result_BA-result_AS);
Err_3=abs(result_ABP-result_AS);
num_1=num_abp;
num_2=num_ba_in;
num_3=num_as;
num_4=num_ba_out;

end