%% Part1 計算像座標真值
dcoord = readtable('UAV_Coordinates_and_Angles.xlsx', 'Sheet', 'Delta_Coordinates', 'VariableNamingRule', 'preserve');
ang = readtable('UAV_Coordinates_and_Angles.xlsx', 'Sheet', 'Angles', 'VariableNamingRule', 'preserve');
%初始已知(units:m)
Xp=184116.126;
Yp=2526216.611;
Zp=330.000;
f=0.05;%(m)
x0=0;
y0=0;

% 設定亂數生成的參數
n = 744;       % 隨機數數量
sigma = 0.02;    % 標準偏差 (單位: mm)
% 生成符合高斯分布 N(0, σ) 的隨機數
random_numbers = sigma * randn(n, 1);

data = cell(25, 5);
data{1, 1} = '片號'; 
data{1, 2} = 'xp(mm)';   
data{1, 3} = 'yp(mm)';  
data{1, 4} = 'xp(mm)(e)';   
data{1, 5} = 'yp(mm)(e)';
a=1;
b=1;

A = zeros(48, 3);
B1 = zeros(48, 3);
L = zeros(48, 1); 
LL=zeros(48, 1);
P = eye(48);
P=P/(sigma*sigma)*1000000;
Start=zeros(24, 3);  
for i=1:8

    omega = double(ang{i, 2}) * (pi / 180);  % 轉換為弧度
    phie = double(ang{i, 3}) * (pi / 180);   % 轉換為弧度
    kappa = double(ang{i, 4}) * (pi / 180);  % 轉換為弧度

    % 計算 m_ij 元素
    m11 = cos(phie) * cos(kappa);
    m12 = sin(omega) * sin(phie) * cos(kappa) + cos(omega) * sin(kappa);
    m13 = -cos(omega) * sin(phie) * cos(kappa) + sin(omega) * sin(kappa);

    m21 = -cos(phie) * sin(kappa);
    m22 = -sin(omega) * sin(phie) * sin(kappa) + cos(omega) * cos(kappa);
    m23 = cos(omega) * sin(phie) * sin(kappa) + sin(omega) * cos(kappa);
    
    m31 = sin(phie);
    m32 = -sin(omega) * cos(phie);
    m33 = cos(omega) * cos(phie);

    for j=1:3
        f=0.05;
        a=a+1;
        %1
        row=(i-1)*3+j;
        delta_x = double(dcoord{row, 3});   
        delta_y = double(dcoord{row, 4});   
        delta_z = double(dcoord{row, 5});
        x_pic=x0-f*(((m11*delta_x)+(m12*delta_y)+(m13*delta_z))/ ...
            ((m31*delta_x)+(m32*delta_y)+(m33*delta_z)));
        y_pic=y0-f*(((m21*delta_x)+(m22*delta_y)+(m23*delta_z))/ ...
            ((m31*delta_x)+(m32*delta_y)+(m33*delta_z)));
       
        
        ch=char(64+i);
        x_pic_error = x_pic*1000 + random_numbers(15 + (a - 2) * 30); % 第15, 45, 75, ...的隨機數
        y_pic_error = y_pic*1000 + random_numbers(30 + (a - 2) * 30); % 第30, 60, 90, ...的隨機數
        L(2*b-1,:)=[x_pic_error];
        L(2*b,:)=[y_pic_error];
        data{a, 1} = sprintf('%s%d', ch, j); 
        data{a, 2} = sprintf('%.3f', x_pic * 1000); % 單位轉為 mm, 保留三位小數
        data{a, 3} = sprintf('%.3f', y_pic * 1000); % 單位轉為 mm, 保留三位小數
        data{a, 4} = sprintf('%.3f', x_pic_error);  
        data{a, 5} = sprintf('%.3f', y_pic_error); 
        b=b+1;
    end
    %未知數近似值
    %像座標
    f=50;
    xpl=L(6*i-5);
    ypl=L(6*i-4);
    xpc=L(6*i-3);
    ypc=L(6*i-2);
    xpr=L(6*i-1);
    ypr=L(6*i);
    %地面座標
    XLl=184116.126-dcoord{(i-1)*3+1, 3};
    YLl=2526216.114-dcoord{(i-1)*3+1, 4};
    ZLl=330-dcoord{(i-1)*3+1, 5};
    XLc=184116.126-dcoord{(i-1)*3+2, 3};
    YLc=2526216.114-dcoord{(i-1)*3+2, 4};
    ZLc=330-dcoord{(i-1)*3+1, 5};
    XLr=184116.126-dcoord{(i-1)*3+3, 3};
    YLr=2526216.114-dcoord{(i-1)*3+3, 4};
    ZLr=330-dcoord{(i-1)*3+1, 5};
    %各參數
    H =330-dcoord{3*i,5};
    B=sqrt(((XLl-XLc)*(XLl-XLc))+((YLl-YLc)*(YLl-YLc)));
    Pp1=xpl-xpc;
    Pp2=xpc-xpr;
    Pp3=xpl-xpr;
    Xp01=B-(xpl/Pp1);
    %(((XLc-XLl)/B)*xpc*10^(-3))-(((YLc-YLl)/B)*ypc*10^(-3))+XLl;
    Yp01=(((XLc-XLl)/B)*ypc*10^(-3))+(((YLc-YLl)/B)*xpc*10^(-3))+YLl;
    Zp01=H-((B*f)/Pp1);
    Start(3*(i-1)+1,:)=[Xp01,Yp01,Zp01];
    Xp02=B-(xpc/Pp2);
    %(((XLr-XLc)/B)*xpr*10^(-3))-(((YLr-YLc)/B)*ypr*10^(-3))+XLc;
    Yp02=(((XLr-XLc)/B)*ypr*10^(-3))+(((YLr-YLc)/B)*xpr*10^(-3))+YLc;
    Zp02=H-((B*f)/Pp2);
    Start(3*(i-1)+2,:)=[Xp02,Yp02,Zp02];
    Xp03=B*2-(xpr/Pp3);
    %=(((XLr-XLl)/(B*2))*xpr*10^(-3))-(((YLr-YLl)/(B*2))*ypr*10^(-3))+XLl;
    Yp03=(((XLr-XLl)/(B*2))*ypr*10^(-3))+(((YLr-YLl)/(B*2))*xpr*10^(-3))+YLl;
    Zp03=H-((B*2*f)/Pp3);
    Start(3*(i-1)+3,:)=[Xp03,Yp03,Zp03];

    q_1=m31*(Xp01-XLl)+m32*(Yp01-YLl)+m33*(Zp01-ZLl);
    r_1=m11*(Xp01-XLl)+m12*(Yp01-YLl)+m13*(Zp01-ZLl);
    s_1=m21*(Xp01-XLl)+m22*(Yp01-YLl)+m23*(Zp01-ZLl);
    b14_1=(f/(q_1*q_1))*(r_1*m31-q_1*m11);
    b15_1=(f/(q_1*q_1))*(r_1*m32-q_1*m12);
    b16_1=(f/(q_1*q_1))*(r_1*m33-q_1*m13);
    b24_1=(f/(q_1*q_1))*(s_1*m31-q_1*m21);
    b25_1=(f/(q_1*q_1))*(s_1*m32-q_1*m22);
    b26_1=(f/(q_1*q_1))*(s_1*m33-q_1*m23);
    
    B1(6*(i-1)+1,1) = b14_1;
    B1(6*(i-1)+1,2) = b15_1;
    B1(6*(i-1)+1,3) = b16_1;
    B1(6*(i-1)+2,1) = b24_1;
    B1(6*(i-1)+2,2) = b25_1;
    B1(6*(i-1)+2,3) = b26_1;
    j1=xpl+f*(r_1/q_1);
    k1=ypl+f*(s_1/q_1);
    LL(6*(i-1)+1, 1)=j1;
    LL(6*(i-1)+2, 1)=k1;

    q_2=m31*(Xp02-XLc)+m32*(Yp02-YLc)+m33*(Zp02-ZLc);
    r_2=m11*(Xp02-XLc)+m12*(Yp02-YLc)+m13*(Zp02-ZLc);
    s_2=m21*(Xp02-XLc)+m22*(Yp02-YLc)+m23*(Zp02-ZLc);
    b14_2=(f/(q_2*q_2))*(r_2*m31-q_2*m11);
    b15_2=(f/(q_2*q_2))*(r_2*m32-q_2*m12);
    b16_2=(f/(q_2*q_2))*(r_2*m33-q_2*m13);
    b24_2=(f/(q_2*q_2))*(s_2*m31-q_2*m21);
    b25_2=(f/(q_2*q_2))*(s_2*m32-q_2*m22);
    b26_2=(f/(q_2*q_2))*(s_2*m33-q_2*m23);
    B1(6*(i-1)+3, :) = [b14_2, b15_2, b16_2];
    B1(6*(i-1)+4, :) = [b24_2, b25_2, b26_2];
    j2=xpc+f*(r_2/q_2);
    k2=ypc+f*(s_2/q_2);
    LL(6*(i-1)+3, 1)=j2;
    LL(6*(i-1)+4, 1)=k2;

    q_3=m31*(Xp03-XLr)+m32*(Yp03-YLr)+m33*(Zp03-ZLr);
    r_3=m11*(Xp03-XLr)+m12*(Yp03-YLr)+m13*(Zp03-ZLr);
    s_3=m21*(Xp03-XLr)+m22*(Yp03-YLr)+m23*(Zp03-ZLr);
    b14_3=(f/(q_3*q_3))*(r_3*m31-q_3*m11);
    b15_3=(f/(q_3*q_3))*(r_3*m32-q_3*m12);
    b16_3=(f/(q_3*q_3))*(r_3*m33-q_3*m13);
    b24_3=(f/(q_3*q_3))*(s_3*m31-q_3*m21);
    b25_3=(f/(q_3*q_3))*(s_3*m32-q_3*m22);
    b26_3=(f/(q_3*q_3))*(s_3*m33-q_3*m23);
    B1(6*(i-1)+5, :) = [b14_3, b15_3, b16_3];
    B1(6*(i-1)+6, :) = [b24_3, b25_3, b26_3];
    j3=xpr+f*(r_3/q_3);
    k3=ypr+f*(s_3/q_3);
    LL(6*(i-1)+5, 1)=j3;
    LL(6*(i-1)+6, 1)=k3;

   
end
fprintf('Start matrix with 6 decimal precision:\n');
for k = 1:size(Start, 1)
    fprintf('%.6f %.6f %.6f\n', Start(k, 1), Start(k, 2), Start(k, 3));
end
%法方程式、未知數近似
N=transpose(B1)*P*B1;
U=transpose(B1)*P*LL;

G=chol(N, 'upper');  % 計算上三角矩陣 G，使得 N = G' * G
Y=G' \ U;
delta=G\Y;
disp(delta);
%% Part2 繪製744個變數的圖
% 生成 744 個隨機數
% 將隨機數標準化為以「幾個標準差」為單位
standardized_numbers = random_numbers / sigma;

% 繪製分布直方圖
figure;
histogram(standardized_numbers, 'Normalization', 'pdf');
hold on;

% 加入理論的正態分布曲線，以「幾個標準差」為單位
x_values = linspace(min(standardized_numbers), max(standardized_numbers), 100); % 橫軸範圍（以標準差為單位）
pdf_values = (1 / sqrt(2 * pi)) * exp(-(x_values .^ 2) / 2);  % N(0, 1) 分布

% 繪製理論正態分布曲線
plot(x_values, pdf_values, 'r-', 'LineWidth', 1.5);

% 設定圖表標題和標籤
title('隨機亂數的分布直方圖 (標準化單位：標準差)');
xlabel('偏差 (標準差的倍數)');
ylabel('頻率');

% 顯示圖表
hold off;

%% Part3
% 讀取 Excel 檔案中的數據
dcoord = readtable('UAV_Coordinates_and_Angles.xlsx', 'Sheet', 'Delta_Coordinates', 'VariableNamingRule', 'preserve');

% 提取 X, Y, Z 座標
X = dcoord{:, 3}; % 第三列為 X 座標
Y = dcoord{:, 4}; % 第四列為 Y 座標
Z = dcoord{:, 5}; % 第五列為 Z 座標

% 繪製三維線條
figure;
hold on;
grid on;

% 確保總點數是三的倍數，否則忽略最後剩餘的點
num_points = floor(length(X) / 3) * 3;

% 每三個點繪製一條線
for i = 1:3:num_points
    % 提取三個點的座標
    X_line = X(i:i+2);
    Y_line = Y(i:i+2);
    Z_line = Z(i:i+2);
    
    % 繪製線條
    plot3(X_line, Y_line, Z_line, '-o', 'LineWidth', 1.5);
end

xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
title('3D Lines from UAV Data');
hold off;
