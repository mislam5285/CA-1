clc;
clear all;
N = 100;
system = 100;
step = 50;
%��ͼ�������
Num1 = zeros(1,step);
Num2 = zeros(1,step);
Num3 = zeros(1,step);
Num4 = zeros(1,step);
for sys = 1:system
   
C=ones(N,N);
ID = zeros(N*N,12);
IDstore = zeros(N*N,12);%�������
%p,m,n,k,l:���ñ���;
%��ʼ��ID����
%���߻���
ID(:,1) = [1:(N*N)]';%ID
for i = 1:N
    ID(i*N-N+1:i*N,2) = i*ones(N,1);
    ID(i*N-N+1:i*N,3) = 1:N; 
end%��γ
for i=1:N*N
    %����
    p = rand(1);
    if p <0.5
        ID(i,4) = 0;%0-14
    end
    if p >= 0.5 && p <= 0.8
        ID(i,4) = 1;%15-44
    end
    if p > 0.8
        ID(i,4) = 2;%45+
    end
    %�Ա�
    p = rand(1);
    if p <=  0.5
        ID(i,5) = 1;% ����
    end 
    ID(i,6) = 0;
    %������
    ID(i,7) =  0.6;
    if ID(i,4) == 0
        ID(i,7) =  0.4;
    end
    %%%%%%%%%%%%%%%%%��Ⱦ��δ���壬��Ҫ��ÿһ�����з������壡
    ID(i,9) =  0;
    ID(i,11) = 0;
    ID(i,12) = 0;
    %��������
    p = rand(1);
    if p<0.1
        ID(i,10) = 1; 
    end
end
%��ʼ��C����
for i=1:N
    for j=1:N
        C(i,j) = N*(i-1)+j;
    end
end
%��ʼ�ǽ�����Ⱥ��Ϣ
%��ʼ��Ǳ����
pp = randi(3)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 1;%Ǳ��
        ID(p,8) = 1;%%%%%%%%%%%%%%%����Ҫ�Ľ���������������
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = 1;
    end
end
%�������
pp = randi(15)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 2;%����
        ID(p,8) = 1;%%%%%%%%%%%%%%%����Ҫ�Ľ���������������
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = randi(10)+7;
    end
end
%���������
pp = randi(7)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 3;%����
        ID(p,8) = 0.9;%%%%%%%%%%%%%%%����Ҫ�Ľ���������������
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,11) = randi(1);%%%%%%%%���������д����ǣ���������
        ID(p,12) = randi(11)+14;
    end
end
%�������
pp = randi(15)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 4;%����
        ID(p,8) = 0.9;%%%%%%%%%%%%%%%����Ҫ�Ľ���������������
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = randi(11)+14;
    end
end
%��ͼ�������
%Num1 = zeros(1,step);
%Num2 = zeros(1,step);
%Num3 = zeros(1,step);
%������Ȭ���������������;
for v=1:step%����
    %��ͬ���߽�������
    if v > 159
        for i = 1:576
            p = randi(N*N);
            tag = 0;
            while ID(p,6) ~= 1 || ID(p,10) == 1 %Ǳ����
            %while ID(p,6) ~= 2 || ID(p,10) == 1%������
            %while ID(p,6) ~= 3|| ID(p,10) == 1%����
                ID(p,12) = 0;
                ID(p,10) = 1;
                ID(p,9) = 0;
                ID(p,6) = 0;
                tag = tag + 1;
                if tag > 500
                    break
                end 
            end
        end
    end
for k=1:N
    for l=1:N
        num = C(k,l);%��ȡ�ø���ID����; 
        if ID(num,6) == 0%��������
            %������Χ������κ�Ӱ��
        end
        if ID(num,6) == 1%Ǳ����Ⱥ
            for m = k-1:k+1
                for n = l-1:l+1%�ھ�����(m,n)
                    if m ~= k && n~= l && m <= N && m >= 1 && n <= N && n >= 1
                        if ID(C(m,n),10) == 0%δ�Բ�������
                            p = rand(1);
                            if p < ID(C(k,l),8)*(1-ID(C(m,n),7))%����Ⱦ
                                if ID(C(m,n),6) == 0%�ھӽ���
                                    ID(C(m,n),6) = 1;
                                    ID(C(m,n),9) = 1;
                                end
                                if ID(C(m,n),6) == 1%�ھ�Ǳ����
                                    if ID(C(m,n),11) == 0
                                        ID(C(m,n),12) = ID(C(m,n),12) + 1;%���ٸ�Ⱦ��
                                    end
                                end
                                if ID(C(m,n),6) == 2%�ھ�Ϊ����
                                    %�����������ĸ�Ⱦ�����κ�Ӱ�죻
                                end
                                if ID(C(m,n),6) == 3%�ھӱ�����
                                    %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
                                    %%%%%%%%%�ɸĽ�������������
                                end
                                if ID(C(m,n),6) == 4%�ھ�Ϊ����
                                    %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
                                    %%%%%%%%%�ɸĽ�������������
                                end
                            end
                        end
                    end
                end
            end
        end
        if ID(num,6) == 2%����
             for m = k-1:k+1
                for n = l-1:l+1%�ھ�����(m,n)
                    if m ~= k && n~= l && m <= N && m >= 1 && n <= N && n >= 1
                        if ID(C(m,n),10) == 0%δ�Բ�������
                            p = rand(1);
                            if p < ID(C(k,l),8)*(1-ID(C(m,n),7))%����Ⱦ����
                                if ID(C(m,n),6) == 0%�ھӽ���
                                    ID(C(m,n),6) = 1;
                                    ID(C(m,n),9) = 1;
                                    ID(C(m,n),12) = ID(C(m,n),12) + 1;%���ٸ�Ⱦ��
                                end
                                if ID(C(m,n),6) == 1%�ھ�Ǳ����
                                    if ID(C(m,n),11) == 0%����δ������
                                        ID(C(m,n),12) = ID(C(m,n),12) + 2;%���ٸ�Ⱦ��
                                    end
                                end
                                if ID(C(m,n),6) == 2%�ھ�Ϊ����
                                    ID(C(m,n),12) = ID(C(m,n),12) + 1;%���ٸ�Ⱦ��
                                    %%%%%%%%%�ɸĽ�������������
                                end
                                if ID(C(m,n),6) == 3%�ھӱ�����
                                    %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
                                    %%%%%%%%%�ɸĽ�������������
                                end
                                if ID(C(m,n),6) == 4%�ھ�Ϊ����
                                    %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
                                    %%%%%%%%%�ɸĽ�������������
                                end
                            end
                        end
                    end
                end
             end  
        end
        if ID(num,6) == 3%����
            %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
            %%%%%%%%%�ɸĽ�������������
        end
        if ID(num,6) == 4%����
            %��������ھӵĸ�Ⱦ�����κ�Ӱ�죻
            %%%%%%%%%�ɸĽ�������������
        end  
        %{
        
        %��������
        if v > 159
            if (ID(num,6) == 2 || ID(num,6) == 1) && (ID(num,12) < 20 || ID(num,11) == 1)%����
                if randi(2) == 1%50%���ʽ�������
                    ID(num,12) == 0;
                    ID(num,10) == 1;
                    ID(num,9) == 0;
                    ID(num,6) == 0;
                end
            end
        end
        %}     
    end
end
%���¸�Ⱦʱ�䣬�ж��Ƿ�任״̬;
for k=1:N
    for l=1:N
        %���¸�Ⱦʱ��
        if ID(C(k,l),6) == 1 || ID(C(k,l),6) == 2 || ID(C(k,l),6) == 3
            if ID(C(k,l),11) == 0
                ID(C(k,l),12) = ID(C(k,l),12) + 1;
            end
            if ID(C(k,l),11) == 1
                p = rand(1);
                if p <0.5 %%%%%%%�����д���ȶ������������
                    ID(C(k,1),12) = ID(C(k,l),12) - 1;
                end
                if ID(C(k,l),12) == 0 %δ�����������ָ������Ļ���
                    ID(C(k,l),6) = 0;%�ָ�����״̬
                    ID(C(k,l),9) = 0;
                    ID(C(k,l),10) = 1;%��ת��Բ�������
                end
            end
        end
        %�жϣ�������״̬
        if ID(C(k,l),6) == 1 && ID(C(k,l),12)-7 == randi(3)%Ǳ���ڣ�����8-10��
            ID(C(k,l),6) = 2;
        end
        if ID(C(k,l),6) == 2%������
            p = rand(1);
            if p < 0.5%�����룬���������д����ǣ�������������
                ID(C(k,l),6) = 3;
            end
            if  ID(C(k,l),12)-14 == randi(11)%����15-25��
            if ID(C(k,l),4) == 0%0-14
                p = rand(1);
                if p < 0.6%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.6%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            if ID(C(k,l),4) == 1%15-45
                p = rand(1);
                if p < 0.45%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.45%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            if ID(C(k,l),4) == 2%45+
                p = rand(1);
                if p < 0.8%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.8%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            end
        end 
        if ID(C(k,l),6) == 3%������
            if  ID(C(k,l),12)-14 == randi(11)%����15-25��
            if ID(C(k,l),4) == 0 && ID(C(k,1),11) == 0%0-14
                p = rand(1);
                if p < 0.6%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.6%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            if ID(C(k,l),4) == 1 && ID(C(k,1),11) == 0%15-45
                p = rand(1);
                if p < 0.45%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.45%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            if ID(C(k,l),4) == 2 && ID(C(k,1),11) == 0%45+
                p = rand(1);
                if p < 0.8%������60%
                    ID(C(k,l),6) = 4;%����
                end
                if p >= 0.8%δ������
                    ID(C(k,l),11) = 1;%����õ�����   
                end
            end
            end 
        end
    end
end
%�жϲ�����������Ⱥ�ƶ�������Ԫ��ϵͳ
if fix(v/40) == v/40
    %����������
    p = randi(N*N);
    while ID(p,6) ~= 0 || ID(p,10) == 1
        p = randi(N*N);
    end
    ID(p,6) = 2;%����
    ID(p,8) = 0.9;%%%%%%%%%%%%%%%����Ҫ�Ľ���������������
    ID(p,9) = 1;
    ID(p,10) = 0;
    ID(p,12) = 7;
end
%������Ⱥ�ƶ�������Ԫ��ϵͳ;��Ԫ���ֳ�4����
if fix(v/3) == v/3 
    C1 = reshape(C(1:N/2,1:N/2),N*N/4,1);
    C2 = reshape(C(1:N/2,N/2+1:N),N*N/4,1);
    C3 = reshape(C(N/2+1:N,1:N/2),N*N/4,1);
    C4 = reshape(C(N/2+1:N,N/2+1:N),N*N/4,1); 
    C1 = reshape(C1(randperm(length(C1))),N/2,N/2);
    C2 = reshape(C2(randperm(length(C2))),N/2,N/2);
    C3 = reshape(C3(randperm(length(C3))),N/2,N/2);
    C4 = reshape(C4(randperm(length(C4))),N/2,N/2);
    C = [C1 C2;C3 C4];
end
Patient = C;
for i = 1:N
    for j = 1:N
        Patient(i,j) = ID(C(i,j),6);
        if ID(C(i,j),6) == 4%����
            Num1(1,v) = Num1(1,v) + 1;
        end
        if ID(C(i,j),6) == 3%����
            Num2(1,v) = Num2(1,v) + 1;
        end
        if ID(C(i,j),6) == 2%����
            Num3(1,v) = Num3(1,v) + 1;
        end
        if ID(C(i,j),6) == 1%Ǳ��
            Num4(1,v) = Num4(1,v) + 1;
        end
    end
end
imagesc(Patient);
pause(0.00001);
%��Ⱦ�˲����˼ǵ�������û�м�1���ظ�����1Ŷ������������������   
end
end
v = 1:step;
plot(v,Num1);%����
title('total death');
figure 
plot(v,Num2);%����
title('ge li');
figure 
plot(v,Num3);%����
title('patient');
plot(v,Num4);%Ǳ��
title('qian fu');

