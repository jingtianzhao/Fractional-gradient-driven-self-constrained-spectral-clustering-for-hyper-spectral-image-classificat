clc;
clear all;
close all;
load indian_pines.mat
load indian_pines_gt.mat

% figure(1);
% [rgb] = func_hyperImshow(indian_pinespart,[10,100,180]);
Hs=indian_pines;

[m,n,r]=size(Hs);
L=m*n;

% R2=zeros(m,n,r);

for i=1:m
    for j=1:n
%         if i==1||i==m||j==1||j==n
        if i==1||i==m||j==1||j==n||i==2||i==m-1||j==2||j==n-1
            HS(i,j,:)=Hs(i,j,:);
        else
            HS(i,j,:)=(Hs(i,j,:)+Hs(i,j-1,:)+Hs(i,j+1,:)+Hs(i,j-2,:)+Hs(i,j+2,:)...
                +Hs(i-1,j,:)+Hs(i-1,j-1,:)+Hs(i-1,j+1,:)+Hs(i-1,j-2,:)+Hs(i-1,j+2,:)...
                +Hs(i+1,j,:)+Hs(i+1,j-1,:)+Hs(i+1,j+1,:)+Hs(i+1,j-2,:)+Hs(i+1,j+2,:)...
                +Hs(i-2,j,:)+Hs(i-2,j-1,:)+Hs(i-2,j+1,:)+Hs(i-2,j-2,:)+Hs(i-2,j+2,:)...
                +Hs(i+2,j,:)+Hs(i+2,j-1,:)+Hs(i+2,j+1,:)+Hs(i+2,j-2,:)+Hs(i+2,j+2,:))/25;
        end
    end
end
        
data=reshape(HS,L,r);
c=100;
alpha=0.1;
beta=0.01;
O=20;
T=20;
tt=clock;

Z=pdist(data);
W=squareform(Z);
 W = W.*W;   
    W = -W/(2*sigma*sigma);
    S = full(spfun(@exp, W));
    D = full(sparse(1:m, 1:m, sum(S))); 
    A =(D^(-1/2) * S * D^(-1/2));

tp2=ones(L,1);
z=diag(tp2,0);
% 初始v
 
for i=1:L
    xxs2=data(i,:);
    xxs1=ews(xxs2);
    xxs(i)=xxs1;
end
B=sort(xxs(:));
B2=B(L-c+1:L,1);
for i=1:c
[D(i),F(i)]=find(B==B2(i));
end
    v=z(:,D);
R=W;
o=0;
yt=0.5;
    h=v;
while o<O   
    % H
    t=0;
    while t<T
        H=(A+2*alpha*(R+lmd)*z)*h/(1+beta)+(beta*v)/(1+beta);
        t=t+1;
        h=H;
    end

    Z=sqrt((H*H').^2./(ones(L,1)*sum((H*H').^2)))/((R-lmd));
    ZZ=z+exp(-o)*Z;
    lmd=yx(ZZ,Z);
    V4=V3;
    V3=V2;
    V2=V1;
    V1=v;
    V=sqrt(H.^2./(ones(L,1)*sum(H.^2)));
    VV=cv(yita,V1,V2,V3,V4,V);
    % re
    z=ZZ;
    v=VV;
    o=o+1;
end 
K=16;
load U2.mat
[idx,ctrs] = kmeans(H,K);%ctrs为中心点，idx为类别数(列矩阵)。
%%
% idx = kmeans(H,4);
[m2,n2]=size(indian_pines_gt);
s=reshape(indian_pines_gt,m2*n2,1);
% idx(find(indian_pines_gtpart==0))=0;
idx = bestMap(s,idx);
%%
c0=find(idx==0);
c1 = find(idx==1);
c2 = find(idx==2);
c3 = find(idx==4);
c4 = find(idx==3);
c5 = find(idx==5);
c6 = find(idx==6);
c7 = find(idx==7);
c8 = find(idx==8);
c9 = find(idx==9);
c10 = find(idx==10);
c11 = find(idx==11);
c12 = find(idx==12);
c13 = find(idx==13);
c14 = find(idx==14);
c15 = find(idx==15);
c16 = find(idx==16);
Z1=zeros(c+1,3);
c0=c0';
[u0,t0]=size(c0);
for s0=1:t0
    Z1(c0(s0),1)= 0;
    Z1(c0(s0),2)= 0;
    Z1(c0(s0),3)= 120;
end
c1=c1';
[u1,t1]=size(c1);
for s1=1:t1
    Z1(c1(s1),1)= 255;
    Z1(c1(s1),2)= 222;
    Z1(c1(s1),3)= 173;
end
c2=c2';
[u2,t2]=size(c2);
for s2=1:t2
    Z1(c2(s2),1)=230;
    Z1(c2(s2),2)=230;
    Z1(c2(s2),3)=250;
end
c3=c3';
[u3,t3]=size(c3);
for s3=1:t3
    Z1(c3(s3),1)=105;
    Z1(c3(s3),2)=105;
    Z1(c3(s3),3)=105;
end
c4=c4';
[u4,t4]=size(c4);
for s4=1:t4
    Z1(c4(s4),1)=0;
    Z1(c4(s4),2)=0;
    Z1(c4(s4),3)=128;
end
c5=c5';
[u5,t5]=size(c5);
for s5=1:t5
    Z1(c5(s5),1)=135;
    Z1(c5(s5),2)=206;
    Z1(c5(s5),3)=235;
end
c6=c6';
[u6,t6]=size(c6);
for s6=1:t6
    Z1(c6(s6),1)=0;
    Z1(c6(s6),2)=255;
    Z1(c6(s6),3)=255;
end
c7=c7';
[u7,t7]=size(c7);
for s7=1:t7
    Z1(c7(s7),1)=0;
    Z1(c7(s7),2)=100;
    Z1(c7(s7),3)=0;
end
c8=c8';
[u8,t8]=size(c8);
for s8=1:t8
    Z1(c8(s8),1)=0;
    Z1(c8(s8),2)=255;
    Z1(c8(s8),3)=127;
end
c9=c9';
[u9,t9]=size(c9);
for s9=1:t9
    Z1(c9(s9),1)=154;
    Z1(c9(s9),2)=205;
    Z1(c9(s9),3)=50;
end
c10=c10';
[u10,t10]=size(c10);
for s10=1:t10
    Z1(c10(s10),1)=255;
    Z1(c10(s10),2)=215;
    Z1(c10(s10),3)=0;
end
c11=c11';
[u11,t11]=size(c11);
for s11=1:t11
    Z1(c11(s11),1)=205;
    Z1(c11(s11),2)=92;
    Z1(c11(s11),3)=92;
end
c12=c12';
[u12,t12]=size(c12);
for s12=1:t12
    Z1(c12(s12),1)=139;
    Z1(c12(s12),2)=69;
    Z1(c12(s12),3)=19;
end
c13=c13';
[u13,t13]=size(c13);
for s13=1:t13
    Z1(c13(s13),1)=255;
    Z1(c13(s13),2)=140;
    Z1(c13(s13),3)=0;
end
c14=c14';
[u14,t14]=size(c14);
for s14=1:t14
    Z1(c14(s14),1)=255;
    Z1(c14(s14),2)=20;
    Z1(c14(s14),3)=147;
end
c15=c15';
[u15,t15]=size(c15);
for s15=1:t15
    Z1(c15(s15),1)=46;
    Z1(c15(s15),2)=139;
    Z1(c15(s15),3)=87;
end
c16=c16';
[u16,t16]=size(c16);
for s16=1:t16
    Z1(c16(s16),1)=0;
    Z1(c16(s16),2)=139;
    Z1(c16(s16),3)=139;
end
Y=reshape(Z1,m,n,3);
y=Y/255;
figure(2)
imshow(y)
%%
% OA 
% idx(find(s==0))=[];
% s(find(s==0))=[];
t=0;
L=length(idx);
for i=1:L
    if idx(i)==s(i)
        t=t+1;
    end
end
  OA=t/(m*n); 
  
%%
disp(['运行时间:  ',num2str(etime(clock,tt))]);