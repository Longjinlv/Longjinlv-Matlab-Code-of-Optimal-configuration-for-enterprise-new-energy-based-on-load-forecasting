clear all
clc;clear

global A fuhe year month day daylenth H0 H1 H2 H3 p2 p3 p4 p z;
path = 'C:\Users\陶陶\Desktop\临时\excelfile\新建文件夹\excelfile\';
namelist = dir([path,'*.xlsx']);
% namelist=sort(namelist.name);
L = length(namelist);
byq=xlsread('中压台帐.xls');
AA=xlsread('编号对应.xls');
% maxfileno=[];
% for i=461:length(namelist)
%     filename{i} = [path,namelist(i).name];
%     D= xlsread(filename{i});%读取excel文件
%     indx=find(abs( D(1,7)-byq(:,1))<0.001);
%     if length(indx)>0
%         maxfileno(i)=indx;
%     end
%     i
% end
%     
% [maxfi,I]=sort(maxfileno);
% indx=find(maxfi>0);
% AA=[maxfi(indx)',I(indx)'];


%% 


AllA = cell(1,L);
js=1;
for i = AA(80,2)
    filename{i} = [path,namelist(i).name];
    [D,DT] = xlsread(filename{i});%读取excel文件
    A=D(:,[7:11,2,13]);
    no=1; y=2;m=3;d=4;t=5;  ap=6; rp=7; % 所在的列位置
    LL=length(A);
    BB=[];
    in=0;
    for year=2022:2023
        for month=1:12
            for day=1:31
                index=find( (abs(A(:,y)-year)<10*eps).*(abs(A(:,m)-month)<10*eps).*(abs(A(:,d)-day)<10*eps));
                for i=1:length(index)
                    if i==1&& index(i)==1
                        k=1;
                    else
                        k=round( (A(index(i),t)-A(index(i)-1,t))*24*60/15);
                    end
                    if k==1
                        in=in+1;
                        BB(in,:)=A(index(i),:);
                    else
                        for j=1:k
                            in=in+1;
                            BB(in,:)=A(index(i),:);
                            BB(in,t)=BB(in-1,t)+15/(24*6);
                        end
                    end
                end
            end
        end
    end
    indx=find( abs(A(1,1)-byq(:,1))<0.1);
    byqlr=byq(indx,2);
end

for i=1:length(BB)
    time=mod((BB(i,t)-24)*24,24);
    if abs(time)<0.0001
        BB(i,7)=24;
    else
        BB(i,7)=time;
    end
    d=BB(i,7);
    
    if d>19+0.01 && d<=21+0.01
        BB(i,8)=1.2857;
    elseif  (d>8 +0.01&& d<=11+0.01) || (d>13+0.01 && d<=19+0.01) || (d>21+0.01 && d<=22+0.01)
        BB(i,8)=1.0739;
    else
        BB(i,8)=0.3434;
    end
end
C=A;
L=length(BB);
ind=find(isnan(BB(:,ap)));
for i=1:length(ind)
    BB(ind(i),ap)=BB(ind(i)-1,ap);
end
BB(L+1,:)=BB(L,:);
BB(L+1,y)=BB(L,y)+1;
A=BB;
%%

H0=500;  % 储能电量，放电最大功率
H1=H0/2; % 储能充电功率
H2=byqlr*1.2;  % 变压器负荷
% tH2=290;
p2=48; % 尖峰电费，基础电费，单价
p3=1500; % 储能成本，每度
p4=100000; % 扩容成本，每100kva
H3=H0/100; % 扩容
ky=1;
z=3 % 年限
dtchushi=H0;
dtchushi00=0.95*H0;
year=2022; month=1;
day=1;
index=[];
no=1; y=2;m=3;d=4;t=5;  ap=6; rp=7; % 所在的列位置
index=find( (  abs(A(:,y)-year)<10*eps ).*( abs(A(:,m)-month)<10*eps).*(abs(A(:,d)-day)<10*eps))+ky;
daylenth=length(index);
fuhe=abs(A(index,ap));
p=A(index,8);  % 电价
L=daylenth;
Agonglv=youhua3(fuhe,dtchushi,100);
dianshouyi1=p'*(Agonglv(:,1)-Agonglv(:,2))/4;

n2=max(fuhe);
mbiao= dianshouyi1(1)+H0*p3/(365*z)*(L/4/24) +max(fuhe+Agonglv(:,1)-Agonglv(:,2))*p2 /30*(L/4/24) -n2*p2 /30*(L/4/24)   -p4*H3/(365*z)*(L/4/24) ;
fuhe1=fuhe;
kb=0.90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tH2=max(fuhe)*kb;
% dtchushi=sum(Agonglv(:,1)-Agonglv(:,2))/4+dtchushi;

%
%
figure
inc=0;
for month=1:12
    dianshouyi=[];
    for day=1:31
        index=[];
        no=1; y=2;m=3;d=4;t=5;  ap=6; rp=7; % 所在的列位置
        index=find( (  abs(A(:,y)-year)<10*eps ).*( abs(A(:,m)-month)<10*eps).*(abs(A(:,d)-day)<10*eps))+ky;
        % index=find( (  abs(A(:,2)-year)<10*eps ).*( abs(A(:,3)-month)<10*eps));
        daylenth=length(index);
        fuhe=abs(A(index,ap));
        hold on
        if daylenth==96
        plot([1:96]/4,fuhe,'.','Color',[rand(),rand(),rand()])
        end
    end
end


figure
inc=0;
for month=1:12
    dianshouyi=[];
    for day=1:31
%         if day==2
%             tH2=max(fuhe)*kb;
%         end
        index=[];
        no=1; y=2;m=3;d=4;t=5;  ap=6; rp=7; % 所在的列位置
        index=find( (  abs(A(:,y)-year)<10*eps ).*( abs(A(:,m)-month)<10*eps).*(abs(A(:,d)-day)<10*eps))+ky;
        % index=find( (  abs(A(:,2)-year)<10*eps ).*( abs(A(:,3)-month)<10*eps));
        daylenth=length(index);
        fuhe=abs(A(index,ap));
        p=A(index,8);  % 电价
        
        if abs(daylenth-96)<0.1
            
            for s=1:96
                if fuhe(s)+Agonglv(s,1)-Agonglv(s,2)>tH2
                    
                    dtchushi2=sum(Agonglv(1:s-1,1)-Agonglv(1:s-1,2))/4+dtchushi;
                    Agonglv(s,1)=0;
                    Agonglv(s,2)=min(fuhe(s)-tH2,dtchushi2);
                    dtchushi2=sum(Agonglv(1:s,1)-Agonglv(1:s,2))/4+dtchushi;

                    mfuhe=[fuhe1(s+1:end)];
                    Agonglv2=youhua3(mfuhe,dtchushi2,dtchushi00);
                    Agonglv=[Agonglv(1:s,:);Agonglv2];
                elseif fuhe(s)+Agonglv(s,1)-Agonglv(s,2)<0
                    Agonglv(s,2)=fuhe(s)+Agonglv(s,1);
                    if s<96
                        dtchushi2=sum(Agonglv(1:s,1)-Agonglv(1:s,2))/4+dtchushi;
                        mfuhe=[fuhe1(s+1:end)];
                        Agonglv2=youhua3(mfuhe,dtchushi2,dtchushi00);
                        Agonglv=[Agonglv(1:s,:);Agonglv2];
                    end
                end
            end
            
            hold on
            plot([1:96]/4,fuhe+Agonglv(:,1)-Agonglv(:,2),'.','Color',[rand(),rand(),rand()])
            
            maxa(day)=max(fuhe);
            maxgong(day)=max(fuhe+Agonglv(:,1)-Agonglv(:,2));
            dianshouyi(day)=p'*(Agonglv(:,1)-Agonglv(:,2))/4;
            dianfei1(day)=p'*fuhe/4;
            dianfei2(day)=p'*(fuhe+Agonglv(:,1)-Agonglv(:,2))/4;
            dtchushi=sum(Agonglv(:,1)-Agonglv(:,2))/4+dtchushi;
            
            
            Agonglv=youhua3(fuhe,dtchushi,dtchushi00);
            fuhe1=fuhe;
            tH2=max(tH2,max(fuhe)*kb);
%              tH2=byqlr;
            
        end
    end
    fval(month)=sum(dianshouyi)+ H0*p3/(12*z)+max(maxgong)*p2  -max(maxa)*p2    -p4*H3/(12*z);
    
    last(month,:)=[sum(dianfei1),max(maxa)*p2,sum(dianfei2),sum(dianfei1)-sum(dianfei2),max(maxgong)*p2,-max(maxgong)*p2+max(maxa)*p2];
    basedianjia(month)=max(maxgong)* p2  -max(maxa)*p2 ;
    shiyongdianjia(month)=sum(dianshouyi);
end

msgbox(['last = ',num2str([sum(fval),sum(basedianjia),sum(shiyongdianjia)])])


