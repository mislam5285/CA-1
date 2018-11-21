clc;
clear all;
N = 100;
system = 100;
step = 50;
%画图所需参数
Num1 = zeros(1,step);
Num2 = zeros(1,step);
Num3 = zeros(1,step);
Num4 = zeros(1,step);
for sys = 1:system
   
C=ones(N,N);
ID = zeros(N*N,12);
IDstore = zeros(N*N,12);%缓存矩阵
%p,m,n,k,l:常用变量;
%初始化ID矩阵
%免疫患者
ID(:,1) = [1:(N*N)]';%ID
for i = 1:N
    ID(i*N-N+1:i*N,2) = i*ones(N,1);
    ID(i*N-N+1:i*N,3) = 1:N; 
end%经纬
for i=1:N*N
    %年龄
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
    %性别
    p = rand(1);
    if p <=  0.5
        ID(i,5) = 1;% 男性
    end 
    ID(i,6) = 0;
    %免疫力
    ID(i,7) =  0.6;
    if ID(i,4) == 0
        ID(i,7) =  0.4;
    end
    %%%%%%%%%%%%%%%%%传染力未定义，需要在每一步长中反复定义！
    ID(i,9) =  0;
    ID(i,11) = 0;
    ID(i,12) = 0;
    %天生免疫
    p = rand(1);
    if p<0.1
        ID(i,10) = 1; 
    end
end
%初始化C矩阵
for i=1:N
    for j=1:N
        C(i,j) = N*(i-1)+j;
    end
end
%初始非健康人群信息
%初始化潜伏者
pp = randi(3)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 1;%潜伏
        ID(p,8) = 1;%%%%%%%%%%%%%%%还需要改进！！！！！！！
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = 1;
    end
end
%随机患者
pp = randi(15)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 2;%患病
        ID(p,8) = 1;%%%%%%%%%%%%%%%还需要改进！！！！！！！
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = randi(10)+7;
    end
end
%随机隔离者
pp = randi(7)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 3;%隔离
        ID(p,8) = 0.9;%%%%%%%%%%%%%%%还需要改进！！！！！！！
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,11) = randi(1);%%%%%%%%被隔离率有待考虑！！！！！
        ID(p,12) = randi(11)+14;
    end
end
%随机死者
pp = randi(15)-1;
if pp >= 1
    for i = 1:pp
        p = randi(N*N);
        ID(p,6) = 4;%死亡
        ID(p,8) = 0.9;%%%%%%%%%%%%%%%还需要改进！！！！！！！
        ID(p,9) = 1;
        ID(p,10) = 0;
        ID(p,12) = randi(11)+14;
    end
end
%画图所需参数
%Num1 = zeros(1,step);
%Num2 = zeros(1,step);
%Num3 = zeros(1,step);
%埃博拉痊愈后具有免疫能力;
for v=1:step%步数
    %不同患者接种疫苗
    if v > 159
        for i = 1:576
            p = randi(N*N);
            tag = 0;
            while ID(p,6) ~= 1 || ID(p,10) == 1 %潜伏者
            %while ID(p,6) ~= 2 || ID(p,10) == 1%患病者
            %while ID(p,6) ~= 3|| ID(p,10) == 1%隔离
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
        num = C(k,l);%获取该个体ID号码; 
        if ID(num,6) == 0%健康个体
            %不对周围人造成任何影响
        end
        if ID(num,6) == 1%潜伏人群
            for m = k-1:k+1
                for n = l-1:l+1%邻居坐标(m,n)
                    if m ~= k && n~= l && m <= N && m >= 1 && n <= N && n >= 1
                        if ID(C(m,n),10) == 0%未对病毒免疫
                            p = rand(1);
                            if p < ID(C(k,l),8)*(1-ID(C(m,n),7))%被传染
                                if ID(C(m,n),6) == 0%邻居健康
                                    ID(C(m,n),6) = 1;
                                    ID(C(m,n),9) = 1;
                                end
                                if ID(C(m,n),6) == 1%邻居潜伏期
                                    if ID(C(m,n),11) == 0
                                        ID(C(m,n),12) = ID(C(m,n),12) + 1;%加速感染；
                                    end
                                end
                                if ID(C(m,n),6) == 2%邻居为患者
                                    %不会加速自身的感染，无任何影响；
                                end
                                if ID(C(m,n),6) == 3%邻居被隔离
                                    %不会加速邻居的感染，无任何影响；
                                    %%%%%%%%%可改进！！！！！！
                                end
                                if ID(C(m,n),6) == 4%邻居为死者
                                    %不会加速邻居的感染，无任何影响；
                                    %%%%%%%%%可改进！！！！！！
                                end
                            end
                        end
                    end
                end
            end
        end
        if ID(num,6) == 2%患者
             for m = k-1:k+1
                for n = l-1:l+1%邻居坐标(m,n)
                    if m ~= k && n~= l && m <= N && m >= 1 && n <= N && n >= 1
                        if ID(C(m,n),10) == 0%未对病毒免疫
                            p = rand(1);
                            if p < ID(C(k,l),8)*(1-ID(C(m,n),7))%被传染病毒
                                if ID(C(m,n),6) == 0%邻居健康
                                    ID(C(m,n),6) = 1;
                                    ID(C(m,n),9) = 1;
                                    ID(C(m,n),12) = ID(C(m,n),12) + 1;%加速感染；
                                end
                                if ID(C(m,n),6) == 1%邻居潜伏期
                                    if ID(C(m,n),11) == 0%病毒未被控制
                                        ID(C(m,n),12) = ID(C(m,n),12) + 2;%加速感染；
                                    end
                                end
                                if ID(C(m,n),6) == 2%邻居为患者
                                    ID(C(m,n),12) = ID(C(m,n),12) + 1;%加速感染；
                                    %%%%%%%%%可改进！！！！！！
                                end
                                if ID(C(m,n),6) == 3%邻居被隔离
                                    %不会加速邻居的感染，无任何影响；
                                    %%%%%%%%%可改进！！！！！！
                                end
                                if ID(C(m,n),6) == 4%邻居为死者
                                    %不会加速邻居的感染，无任何影响；
                                    %%%%%%%%%可改进！！！！！！
                                end
                            end
                        end
                    end
                end
             end  
        end
        if ID(num,6) == 3%隔离
            %不会加速邻居的感染，无任何影响；
            %%%%%%%%%可改进！！！！！！
        end
        if ID(num,6) == 4%死者
            %不会加速邻居的感染，无任何影响；
            %%%%%%%%%可改进！！！！！！
        end  
        %{
        
        %加入疫苗
        if v > 159
            if (ID(num,6) == 2 || ID(num,6) == 1) && (ID(num,12) < 20 || ID(num,11) == 1)%患者
                if randi(2) == 1%50%概率接种疫苗
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
%更新感染时间，判断是否变换状态;
for k=1:N
    for l=1:N
        %更新感染时间
        if ID(C(k,l),6) == 1 || ID(C(k,l),6) == 2 || ID(C(k,l),6) == 3
            if ID(C(k,l),11) == 0
                ID(C(k,l),12) = ID(C(k,l),12) + 1;
            end
            if ID(C(k,l),11) == 1
                p = rand(1);
                if p <0.5 %%%%%%%概率有待商榷！！！！！！
                    ID(C(k,1),12) = ID(C(k,l),12) - 1;
                end
                if ID(C(k,l),12) == 0 %未死亡并慢慢恢复健康的患者
                    ID(C(k,l),6) = 0;%恢复健康状态
                    ID(C(k,l),9) = 0;
                    ID(C(k,l),10) = 1;%好转后对病毒免疫
                end
            end
        end
        %判断，并更新状态
        if ID(C(k,l),6) == 1 && ID(C(k,l),12)-7 == randi(3)%潜伏期，患病8-10天
            ID(C(k,l),6) = 2;
        end
        if ID(C(k,l),6) == 2%发病期
            p = rand(1);
            if p < 0.5%被隔离，被隔离率有待考虑！！！！！！！
                ID(C(k,l),6) = 3;
            end
            if  ID(C(k,l),12)-14 == randi(11)%患病15-25天
            if ID(C(k,l),4) == 0%0-14
                p = rand(1);
                if p < 0.6%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.6%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            if ID(C(k,l),4) == 1%15-45
                p = rand(1);
                if p < 0.45%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.45%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            if ID(C(k,l),4) == 2%45+
                p = rand(1);
                if p < 0.8%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.8%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            end
        end 
        if ID(C(k,l),6) == 3%隔离期
            if  ID(C(k,l),12)-14 == randi(11)%患病15-25天
            if ID(C(k,l),4) == 0 && ID(C(k,1),11) == 0%0-14
                p = rand(1);
                if p < 0.6%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.6%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            if ID(C(k,l),4) == 1 && ID(C(k,1),11) == 0%15-45
                p = rand(1);
                if p < 0.45%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.45%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            if ID(C(k,l),4) == 2 && ID(C(k,1),11) == 0%45+
                p = rand(1);
                if p < 0.8%病死率60%
                    ID(C(k,l),6) = 4;%死者
                end
                if p >= 0.8%未病死者
                    ID(C(k,l),11) = 1;%病情得到控制   
                end
            end
            end 
        end
    end
end
%判断步长，考虑人群移动，打乱元胞系统
if fix(v/40) == v/40
    %患病者移入
    p = randi(N*N);
    while ID(p,6) ~= 0 || ID(p,10) == 1
        p = randi(N*N);
    end
    ID(p,6) = 2;%患病
    ID(p,8) = 0.9;%%%%%%%%%%%%%%%还需要改进！！！！！！！
    ID(p,9) = 1;
    ID(p,10) = 0;
    ID(p,12) = 7;
end
%考虑人群移动，打乱元胞系统;将元胞分成4部分
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
        if ID(C(i,j),6) == 4%死者
            Num1(1,v) = Num1(1,v) + 1;
        end
        if ID(C(i,j),6) == 3%隔离
            Num2(1,v) = Num2(1,v) + 1;
        end
        if ID(C(i,j),6) == 2%患者
            Num3(1,v) = Num3(1,v) + 1;
        end
        if ID(C(i,j),6) == 1%潜伏
            Num4(1,v) = Num4(1,v) + 1;
        end
    end
end
imagesc(Patient);
pause(0.00001);
%感染了病的人记得想想有没有加1或重复加了1哦！！！！！！！！！   
end
end
v = 1:step;
plot(v,Num1);%死者
title('total death');
figure 
plot(v,Num2);%隔离
title('ge li');
figure 
plot(v,Num3);%患者
title('patient');
plot(v,Num4);%潜伏
title('qian fu');

