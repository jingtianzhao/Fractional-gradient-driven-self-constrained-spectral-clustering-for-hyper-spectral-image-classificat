%%信息熵计算方法
function  xxs1=ews(xxs2)
[M,N]=size(xxs2);
temp=zeros(1,256);
for m=1:M
  for n=1:N
    if A(m,n)==0
      i=1;
    else
      i=A(m,n);
    end
  temp(i)=temp(i)+1;
  end
end
temp=temp/(M*N);
result=0;
for i=1:length(temp)
  if temp(i)==0
    result=result;
  else
    result=result-temp(i)*log2(temp(i));
  end
end
xxs1=result;
