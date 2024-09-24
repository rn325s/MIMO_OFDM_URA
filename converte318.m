
function [C0,C1]=converte318(L)
%[C0,C1]=converte(L)
%����LLR�õ���Ӧ��ֵ
%LLR = ln[p(r|b=0)/p(r|b=1)]
C0=zeros(1,length(L));
for i=1:length(L)
    C0(i)=exp(L(i))/(1+exp(L(i)));
    if C0(i)==0
        C0(i)=1e-10;
        C1(i)=1;
    end
end


C1=zeros(1,length(L));
for i=1:length(L)
    C1(i)=1/(1+exp(L(i)));
    if C1(i)==0
        C1(i)=1e-10;
        C0(i)=1;
    end
    
end
end
%��������Ȼ��ת��Ϊ�ж�Ϊ0Ϊ1�ĸ���