function [Agonglv2,fval1]=youhua3(fuhe,dtchushi,dtchushi00)

global  H0 H1 H2 H3 p2 p3 p4 p z;
L=length(fuhe);
% fuhe=mfuhe;

ln=2*L+1;
%1. n>b, xi-yi-n<=-ai 最大功率，也即小于当天的最大需量
A1=zeros(L,ln); b1=zeros(L,1);
d=1;
% [sc,I]=sort(fuhe,'descend' );
% fuhe2=fuhe;
% fuhe2(I(1:d))=fuhe(I(1));
for i=1:L
    A1(i,i)=1;   A1(i,L+i)=-1;   A1(i,2*L+1)=-1;  b1(i)=-fuhe(i);
end


%2.   yi<=@sum(k<i:x(k)-y(k)))+ 储能初始电量;
A21=zeros(L,ln); A22=zeros(L,ln); b2=ones(L,1)*dtchushi*4;
for i=1:L
    A21(i,i+L)=1;
    for j=1:i-1
        A22(i,j)=-1;    A22(i,L+j)=1;
    end
end
A2=A21+A22;

%3.(xi+c)*1/4  + 储能初始电量<H0,  sum(k<i:x(k)-y(k)));  储能电量负荷
A3=zeros(L,ln); b3=ones(L,1)*H0*4 - dtchushi*4;
for i=1:L
    A3(i,i)=1;
    for j=1:i-1
        A3(i,j)=1; A3(i,L+j)=-1;
    end
end
%4 x<H1,充电功率限制
A4=zeros(L,ln); b4=ones(L,1)*H1;
for i=1:L
    A4(i,i)=1;
end
%5 ai+xi-yi<H2，  变压器功率限制
A5=zeros(L,ln); b5=zeros(L,1);
for i=1:L
    A5(i,i)=1;  A5(i,L+i)=-1;  b5(i)=H2-fuhe(i);
end
%6 y<H0，储能放电功率限制
A6=zeros(L,ln); b6=ones(L,1)*H0;
for i=1:L
    A6(i,i+L)=1;
end

%7 ai+xi-yi>0，储能实际放电限制
A7=zeros(L,ln); b7=zeros(L,1);
for i=1:L
    A7(i,i)=-1; A7(i,i+L)=1; b7(i)=fuhe(i);
end

%8 
A8=zeros(1,ln); b8=dtchushi00;
A8(1:L)=-1; A8(L+1:2*L)=+1;

An=[A1;A2;A3;A4;A5;A6;A7;A8];
bn=[b1;b2;b3;b4;b5;b6;b7;b8];



%6 mub
M1=zeros(1,ln);
for j=1:L
    M1(j)=p(j)/4;    M1(L+j)=-p(j)/4;
end
M1(ln)=p2 /30*(L/4/24) *1;
%
lb=zeros(2*L+1,1);
[XX,fval1]= linprog(M1,An,bn,[],[],lb,[]);
Agonglv2=[ XX(1:L,1),XX(L+1:2*L,1)];