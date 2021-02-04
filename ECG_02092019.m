clear all;
clc;
close all;
%***************************************************************
% Get Data & User Inputs
%***************************************************************

Fileloc = 'F:\TNTECH Mamun MS\Summer 18\matlab';    
Filename = input('Enter ECG File Name = ','s'); % Input Filename
% f_bp = input('Enter bp Name = ','s'); % Input bp
% f_ox = input('Enter ox Name = ','s'); % Input ox
% c_bp=str2num(f_bp)
% c_ox=str2num(f_ox)
c_bp=120;
c_ox=96;

Headerfile = strcat(Filename,'.hea'); % Header In TXT format
load(Filename)
% Load Header Data
%***************************************************************
fprintf(1,'\nK> Loading Data from Header File %s ...\n',Headerfile);

signalh = fullfile(Fileloc, Headerfile);
fid1 = fopen(signalh,'r');
z = fgets(fid1);
A = sscanf(z, '%*s %d %d',[1,2]);
nosig = A(1); % Number Of Signals
sfreq = A(2);

clear A;

z = fgetl(fid1);
A = sscanf(z, '%*s %*d %d %d %d %d',[1,4]);
gain = A(2); % Integers Per mV
clear A;


S = sfreq*10;
counter1=0;
counter2=0;
counter3=0;



j = S*0+1:1:S*(0+1);
D = val(1,:);% 2.5.12.14.25 (ST elevation).26 (ST ele).27 (partially ST ).31 ideal ST .05500001m (ST d, hyper T)


dat = length(D);
k = 1:1:dat;
D1 = D(k)/gain;
QQ=D1;
Tm=(j/sfreq);

figure 
plot(Tm,QQ);xlabel('Time');ylabel('Amplitude');hold on
title('Original ECG signal'); 

% %%%Blood pressure 
% 
% prompt = 'What is the number of row of BP signal? '
% BP_S = input(prompt);
% D_BP = val(BP_S,:)
% % DP=min(D_BP)
% % SP=max(D_BP)

% % % heart rate from PPG starts  %%

prompt = 'What is the number of row of PPG signal? '
PPG_S = input(prompt);
D_PPG = val(PPG_S,:);% 2.5.12.14.25 (ST elevation).26 (ST ele).27 (partially ST ).31 ideal ST 


dat_PPG = length(D_PPG);
k_PPG = 1:1:dat_PPG;
D_PPG1 = D_PPG(k_PPG)/gain;
QQ_PPG=D_PPG1+.64;
% QQ_PPG=D_PPG1+.84;
Tm=(j/sfreq);

figure 
plot(Tm,QQ_PPG);xlabel('Time');ylabel('Amplitude');
title('Original PPG signal');

thresh_PPG = 0.25;
%%create time axis
len_PPG = length(QQ_PPG);
tt_PPG = 1/sfreq:1/sfreq:ceil(len_PPG/sfreq);
t_PPG = tt_PPG(1:len_PPG);
max_h_PPG = max (QQ_PPG(round(len_PPG/4):round(3*len_PPG/4)));%

poss_reg_PPG = QQ_PPG>(thresh_PPG*max_h_PPG); %then build an array of segments to look in
left_PPG = find(diff([0 poss_reg_PPG])==1); % remember to zero pad at start
right_PPG = find(diff([poss_reg_PPG 0])==-1); % remember to zero pad at end
for(i=1:length(left_PPG))
[maxval_PPG(i) maxloc_PPG(i)] = max( QQ_PPG(left_PPG(i):right_PPG(i)) );
[minval_PPG(i) minloc_PPG(i)] = min( QQ_PPG(left_PPG(i):right_PPG(i)) );
maxloc_PPG(i) = maxloc_PPG(i)-1+left_PPG(i); % add offset of present location
minloc_PPG(i) = minloc_PPG(i)-1+left_PPG(i); % add offset of present location
end
R_index_PPG = maxloc_PPG;

R_t_PPG = t_PPG(maxloc_PPG);
R_amp_PPG = maxval_PPG;

for j = 2:length (R_t_PPG)
HR_PPG(j)= R_t_PPG(j)-R_t_PPG(j-1);
end
H_R_PPG = round(60/(mean (HR_PPG)))

QQ_PPGInv=1.01*max(QQ_PPG) - QQ_PPG;
[Minima,MinIdx] = findpeaks(QQ_PPGInv);
Minima = QQ_PPG(MinIdx);



%
S_index_PPG=0;
S_amp_PPG=0;
S_t_PPG=0;
K1_PPG=0;
R_PPG_len=0;
IR1_PPG=0;
K1_PPG=mean(diff(R_index_PPG));
R_PPG_len= length (R_index_PPG);
for j = 1:(R_PPG_len-1)
IR1_PPG = R_index_PPG(j);

for i = IR1_PPG:(IR1_PPG+round(K1_PPG/2))
if round((diff(QQ_PPG(i:i+1))/diff(QQ_PPG(i+6:i+7))),0)==3
    S_index_PPG(j)= i+6;
S_amp_PPG(j) = QQ_PPG(1,i+6);
S_t_PPG(j) = t_PPG(1,i+6);

end
end
end
figure 
plot(t_PPG,QQ_PPG);hold on;
plot(S_t_PPG,S_amp_PPG,'*');hold on;
plot(R_t_PPG,R_amp_PPG,'*'); hold on;
plot(t_PPG(MinIdx),Minima,'*');hold on 
% plot(R_t,R_amp,'+r'); 

% % area calculation
% xx=t_PPG(31:115);
% yx=QQ_PPG(31:115);
% yyx = zeros(1,length(xx));
% yyx(:) = 0 ;
% y1x=0:0.01:.9;
% x1x = ones(1,length(y1x));
% x1x(:) = 0.64;
% plot(xx,yx,xx,yyx,'--',x1x,y1x,'--')
% x_area=abs(trapz(xx,yx))
% 
% xx1=t_PPG(72:115);
% yx1=QQ_PPG(72:115);
% yyx1 = zeros(1,length(xx1));
% yyx1(:) = 0 ;
% y1x1=0:0.01:.9;
% x1x1 = ones(1,length(y1x1));
% x1x1(:) = 0.92;
% plot(xx1,yx1,xx1,yyx1,'--',x1x1,y1x1,'--')
% x_area1=abs(trapz(xx1,yx1))




% [Maxima,MaxIdx] = findpeaks(Data);
% DataInv = 1.01*max(Data) - Data;
% [Minima,MinIdx] = findpeaks(DataInv);
% Minima = Data(MinIdx);

% % heart rate with PPG ends here %



% figure 
% 
% plot(Tm(40:280),QQ(40:280));xlabel('Time');ylabel('Amplitude');
% title('Original ECG signal');

%figure 2
D= transpose (D);
windowSize = 5;
filsig = filter (ones(1,windowSize)/windowSize,1,D);

y = medfilt1(filsig,200); % 1st median filter
s1 = y;
clear y;
y = medfilt1(s1,600); % 2nd median filter
D = filsig - y;
QW=D;

% figure 
% 
% plot(Tm,QW);
% title('After filtering and removal of Baseline wander');
% xlabel('Time');ylabel('Amplitude');
% 
% figure 
% plot(Tm(40:280),QW(40:280));
% title('After filtering and removal of Baseline wander');
% xlabel('Time');ylabel('Amplitude');
clear s1;
clear y;

%figure 3
D = transpose (D);
D = cwt (D, 1:4, 'bior2.4');%Performing Continuous Wavelet
%Transform using Biorthogonal
% Wavelet to ECG_1

D = transpose (D);
x = D (:,4);

% %creating figure 3
% panhandle = uipanel('Position'[.25 .1 .67 .67])
% h1 = subplot(4,1,1,'Parent',pandhandle);
% plot(h1,Tm(200:1200),D(200:1200,1));hold on;
% h2 = subplot(4,1,2,'Parent',pandhandle);
% plot(h2,Tm(200:1200),D(200:1200,2));
% h3 = subplot(4,1,3,'Parent',pandhandle);
% plot(h3,Tm(200:1200),D(200:1200,3));
% % h4 = subplot(4,1,4,'Parent',pandhandle);
% plot(h4,Tm(200:1200),D(200:1200,4));
% xlabel('Time');ylabel('Amplitude');

thresh = 0.6;
% create time axis
len = length (x);
tt = 1/sfreq:1/sfreq:ceil(len/sfreq);
t = tt(1:len);
max_h = max (x(round(len/4):round(3*len/4)));%

poss_reg = x>(thresh*max_h); %then build an array of segments to
% look in
left = find(diff([0 poss_reg'])==1); % remember to zero pad at
% start
right = find(diff([poss_reg' 0])==-1); % remember to zero pad at
% end
for(i=1:length(left))
[maxval(i) maxloc(i)] = max( QW(left(i):right(i)) );
[minval(i) minloc(i)] = min( QW(left(i):right(i)) );
maxloc(i) = maxloc(i)-1+left(i); % add offset of present
% location
minloc(i) = minloc(i)-1+left(i); % add offset of present
% location
end
R_index = maxloc;

R_t = t(maxloc);
R_amp = maxval;

% PTT_f=0;
% for i=1:length(R_t)
%     PTT1 = R_t(i)- R_t_PPG(i);
%     PTT_f=PTT_f+PTT1;
% end
% PTT=PTT_f/length(R_t)

%***************************************************************
% Heart Rate Calculation
%***************************************************************
for j = 2:length (R_t)
HR(j)= R_t(j)-R_t(j-1);
end
H_R = round(60/(mean (HR)))



figure 
plot(t_PPG,QQ_PPG);hold on;
plot(S_t_PPG,S_amp_PPG,'*');hold on;
plot(R_t_PPG,R_amp_PPG,'*'); hold on;
plot(t_PPG(MinIdx),Minima,'*');hold on 
% plot(t,QW,'g');hold on;
plot(R_t,R_amp/1000,'+r'); hold on;


fprintf (1,'\nK> Heart rate is %d \n',H_R)

PTT_f=0;
for i=1:min(length(R_t_PPG),length(R_t))
    PTT1(i) = R_t(i)- R_t_PPG(i);
    PTT_f=PTT_f+PTT1(i);
end
PTT=PTT_f/min(length(R_t_PPG),length(R_t));
    
if PTT<0
    
    PTT_final = abs(R_t(4)-R_t(3)+PTT)
else
    PTT_final = PTT
end
% %***************************************************************
% % S-Point Detection
% %***************************************************************
% R_len= length (R_index);
% for j = 1:R_len
% IR1 = R_index(j);
% for i = IR1:IR1+ (round(sfreq*0.2)*(H_R/72))
% if i == length (x)||i==0
% S_index(j)= 1;
% S_amp(j) = QW(1,1);
% S_t(j) = t(1,1);
% break
% end
% if(i-1)==0
%     break
%     
% elseif QW(i,1)< QW(i+1,1) && QW(i,1)< QW(i-1,1)
% S_index(j)= i;
% S_amp(j) = QW(i,1);
% S_t(j) = t(1,i);
% break
% end
% if QW(i,1)== QW(i+1,1) && QW(i,1)< QW(i-1,1)
% S_index(j)= i;
% S_amp(j) = QW(i,1);
% S_t(j) = t(1,i);
% break
% end
% S_index(j) = i;
% S_amp(j)= QW(i,1);
% S_t(j)=t(1,i);
% end
% end
% 
% 
% %slope point detection from PPG 
% 
% 
% 
% %***************************************************************
% % Q-Point Detection
% %***************************************************************
% for j = 1:R_len
% IR1 = R_index(j);
% for i = IR1:-1:IR1- (round(sfreq*0.06 *(H_R/72)))
% if i == 0|i==length (QW)|i-1==0
% Q_index(j)= 1;
% Q_amp(j) = QW(1,1);
% Q_t(j) = t(1,1);
% break
% end
% if QW(i,1)< QW(i+1,1) && QW(i,1)< QW(i-1,1)
% Q_index(j)= i;
% Q_amp(j) = QW(i,1);
% Q_t(j) = t(1,i);
% break
% end
% if QW(i,1)< QW(i+1,1) && QW(i,1)== QW(i-1,1)
% Q_index(j)= i;
% Q_amp(j) = QW(i,1);
% Q_t(j) = t(1,i);
% break
% end
% Q_index(j) = i;
% Q_amp(j)= QW(i,1);
% Q_t(j)=t(1,i);
% end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %***************************************************************
% % J-Point Detection
% %***************************************************************
% S_len = length (S_index);
% for j = 1:S_len
% IS1= S_index(j);
% 
% for i=IS1:1:IS1+(round(sfreq*0.1*(H_R/72)))
% if i==0 || i==1
% J_index(j)=1;
% J_amp(j)= QW(1,1);
% J_t(j)= t(1,1);
% 
% break
% end
% 
% if QW(i,1)>=0
% J_index(j)=i;
% J_amp(j)= QW(i,1);
% J_t(j)= t(1,i);
% 
% break
% end
% 
% 
% 
% if (i-1)==0 || (i-1)<=0 ||(i+1)>=5000
%     break
% elseif QW(i,1)> QW(i+1,1) && QW(i,1)> QW(i-1,1)
% J_index(j)= i;
% J_amp(j) = QW(i,1);
% J_t(j) = t(1,i);
% break
% end
% if QW(i,1)== QW(i+1,1) && QW(i,1)> QW(i-1,1)
% J_index(j)= i;
% J_amp(j) = QW(i,1);
% J_t(j) = t(1,i);
% break
% end
% 
% 
% J_index(j) = i;
% J_amp(j)= QW(i,1);
% J_t(j)=t(1,i);
% end
% % if i==0
% % J_index(j)=1;
% % J_amp(j)= QW(1,1);
% % J_t(j)= t(1,1);
% % end
% 
% 
% end
% 
% 
% %***************************************************************
% % T-Peak Detection
% %***************************************************************
% % J_len = length (J_index);
% % for j = 1: J_len
% % P1 = round(R_index(j)+ (ceil(sfreq*0.6)*(H_R/72)));
% % 
% % P2 = round(J_index (j)+ (ceil(sfreq*0.08) *(H_R/72)));
% % 
% % 
% % if P1> length (x)||P2> length (x)
% % break
% % end
% % if P1 > P2
% % [T_peak(j),T_peak_index(j)] = max(QW(P2:P1));
% % T_peak_index(j) = T_peak_index(j)+ P2;
% % else
% % [T_peak(j),T_peak_index(j)] = max(QW(P1:P2));
% % T_peak_index(j) = T_peak_index(j)+ P1;
% % end
% % T_peak_t (j)= t(round((T_peak_index (j))));
% % end
% 
% 
% 
% 
% 
% % %***************************************************************
% % % T-Peak Detection
% % %***************************************************************
% J_len = length (J_index);
% for j = 1: J_len
% % P1 = round(J_index(j)+ (ceil(sfreq*0.5)*(H_R/72)));
% P1 = round(J_index(j)+ (R_index(2)-J_index(1))/2);
% P2 = round(J_index (j)+ (ceil(sfreq*0.08) *(H_R/72)));
% 
% 
% if P1> length (x)||P2> length (x)
% break
% end
% if P1 > P2
% [T_peak(j),T_peak_index(j)] = max(QW(P2:P1));
% T_peak_index(j) = T_peak_index(j)+ P2;
% else
% [T_peak(j),T_peak_index(j)] = max(QW(P1:P2));
% T_peak_index(j) = T_peak_index(j)+ P1;
% end
% T_peak_t (j)= t(round((T_peak_index (j))));
% end
% %***************************************************************
% % Hyperacute T wave check
% %***************************************************************
% 
% avg_T_peak=mean(T_peak);
% H_T_peak=(avg_T_peak)/(mean(R_amp));
% 
% % if abs(H_T_peak) >.5
% % 
% % fprintf(1,'H_T_peak yes\n')
% % else
% % 
% % fprintf(1,'H_T_peak no\n')
% % end
% 
% %***************************************************************
% % QT prolongation check
% %***************************************************************
% 
% J_len = length (J_index);
% 
% Q_prol = length(J_index);
% 
% for j = 1: (J_len-1)
%     
%   Q_prol(j)= abs(T_peak_t(j)-Q_t(j));
%   
% end
% 
% mean_path_QT=mean(Q_prol);
% 
% % if mean_path_QT >.4
% % fprintf(1,'\nK> QT interval symptom yes \n');
% % else
% % fprintf(1,'\nK> QT interval symptom no \n');
% % end
% 
% 
% %***************************************************************
% % K-Point
% %***************************************************************
% Q_len = length (Q_index);
% for j = 1:Q_len
% IQ1 = Q_index(j);
% foundk = 0;
% for i=IQ1:-1:IQ1- (round(sfreq*0.02*(H_R/72)))
% if i == 0
% K_index(j) = 1;
% K_amp(j)= QW(1,1);
% K_t(j)=t(1,1);
% foundk = 1;
% break
% end
% if QW(i,1)>=0
% K_index(j) = i;
% K_amp(j)= QW(i,1);
% K_t(j)=t(1,i);
% foundk = 1;
% break
% end
% end
% if foundk == 0
% K_index(j)=i;
% K_amp(j)= QW(i,1);
% K_t(j)=t(1,i);
% end
% end
% 
% % %***************************************************************
% % % P-Peak Detection
% % %***************************************************************
% PP_len = length (K_index);
% for j = 1: PP_len
% % PP1 = round(P_index(j)- (ceil(sfreq*0.2)*(H_R/72)));
% PP1 = round(K_index(j)- ((S_index(j)-Q_index(j))*2.5));
% PP2 = round(K_index(j));
% 
% 
% if PP1<0||PP2<0||PP1==0||PP2==0
% break
% 
% elseif PP1 < PP2
% [P_peak(j),P_peak_index(j)] = max(QW(PP2:-1:PP1));
% if (PP2-P_peak_index(j))<=0
%     P_peak_index(j) = K_index(j);
% else
% P_peak_index(j) = PP2-P_peak_index(j);
% end
% % else
% % [P_peak(j),P_peak_index(j)] = max(QW(PP2:1:PP1));
% % P_peak_index(j) = P_peak_index(j)- PP1;
% if P_peak_index(j)<=0
%     P_peak_t(j)=K_index(j);
% else
% P_peak_t(j)= t(round((P_peak_index(j))));
% end
% 
% 
% end
% end
% 
% 
% %***************************************************************
% % P-Point Detection via K+80ms
% %***************************************************************
% K_len = length (K_index);
% for j = 1:K_len
% IK1 = K_index(j);
% for i=IK1:-1:IK1- (round(sfreq*0.03) *(H_R/72))
% if i ==0
% P_index(j) = 1;
% P_amp(j)= QW(1,1);
% P_t(j)=t(1,1);
% break
% end
% P_index(j) = i;
% P_amp(j)= QW(i,1);
% P_t(j)=t(1,i);
% end
% end
% 
% %***************************************************************
% % T-Point (T onset) Detection via T-peak
% %***************************************************************
% T_len = length (T_peak_index);
% for j = 1:T_len
% IT1 = T_peak_index(j);
% for i = IT1-(S_index(2)-Q_index(2)):-1:J_index(j)
%  
% if(i-1)==0 || (i-1)<=0
%     break
% end 
%  
% if QW(i,1)< QW(i+1,1) && QW(i,1)< QW(i-1,1)
% TP_index(j)= i;
% TP_amp(j) = QW(i,1);
% TP_t(j) = t(1,i);
% break
% end    
%  if QW(i,1)< QW(i+1,1) && QW(i,1)== QW(i-1,1)
% TP_index(j)= i;
% TP_amp(j) = QW(i,1);
% TP_t(j) = t(1,i);
% break
%  end  
% 
%   if QW(i,1)== QW(i+1,1) && QW(i,1)== QW(i-1,1)
% TP_index(j)= i;
% TP_amp(j) = QW(i,1);
% TP_t(j) = t(1,i);
% break
% end 
% % if round(QW(i,1),2)== round(J_amp(j),2)
% % TP_index(j)= i;
% % TP_amp(j) = QW(i,1);
% % TP_t(j) = t(1,i);
% % break
% % end
% 
% TP_index(j)=round(i);
% TP_amp(j)=QW(round(i),1);
% TP_t(j)=t(1,round(i));
% 
% end
% % TP_index(j)=J_index(j)+(K_index(j)-P_index(j));
% % TP_amp(j)=QW(round(TP_index(j)),1);
% % TP_t(j)=t(1,round(TP_index(j)));
% end
% 
% 
% %***************************************************************
% % Pathological Q wave 
% %***************************************************************
% Path_len= length (Q_index);
% for j = 1:Path_len
% 
% path_c(j) = Q_amp(j)/R_amp(j);
% 
% end 
% mean_path_c=mean(path_c);
% 
% % if mean_path_c<=-0.25
% % fprintf(1,'\nK> mean path c yes \n');
% % else
% % fprintf(1,'\nK> mean path c no \n');
% % end
% 
% 
% 
% 
% 
% %***************************************************************
% % Calculation of Isoelectric Line
% %***************************************************************
% for j = 1:1:K_len;
% ISO(j) = mean(QW(P_index(j):K_index(j)));
% end
% 
% 
% 
% %***************************************************************
% % ST depression detection using Isoscale line
% %***************************************************************
% 
% 
% ST_dep= length (Q_index);
% for j = 1:ST_dep
% ST_dep1 = Q_index(j);
% if (Q_index(j)+round(sfreq*0.12))<1000
%     
% ST_damp(j) = QW((Q_index(j)+round(sfreq*0.12)),1);
% ST_dt(j) = t(1,(Q_index(j)+round(sfreq*0.12)));
% 
% ST_dvalue(j)=ISO(j)-ST_damp(j);
% else 
%     break
% end
% end 
% mean_ST_dvalue=mean(ST_dvalue);
% 
% % if mean_ST_dvalue >.1
% %     warn = 'warning:ST_Depressions';
% %     fprintf(1,'\nK> ST depression \n');
% % else
% % %     warn = 'No MI';
% %     fprintf(1,'\nK> no ST depression \n');
% % end
% 
% 
% 
% 
% %***************************************************************
% % Calculation of ST-segment
% %***************************************************************
% a = length (J_index);
% b = length (TP_index);
% if a==b;
% for j = 1:1:J_len;
% ST(j) = mean(QW(J_index(j):TP_index(j)));
% end
% end
% if a>b
% for j = 1:1:b;
% ST(j) = mean(QW(J_index(j):TP_index(j)));
% end
% end
% if a<b
% for j = 1:1:a;
% ST(j) = mean(QW(J_index(j):TP_index(j)));
% end
% end
% %***************************************************************
% % Comparison of ISO and ST
% %***************************************************************
% a = length (ISO);
% b = length (ST);
% if a==b
% for j = 1:a
% counter1=counter1+1;
% if (ISO(j))>= (ST(j)+0.05) || ISO(j)>=ST(j)-0.05||ISO(j)==ST(j)
% 
% counter2=counter2+1;
% else
% counter3=counter3+1;
% end
% end
% 
% end
% if a<b
% for j=1:a
% counter1=counter1+1;
% if ISO(j)>=ST(j)+0.05 || ISO(j)>=ST(j)-0.05||ISO(j)==ST(j)
% counter2=counter2+1;
% else
% counter3=counter3+1;
% end
% end
% end
% if a>b
% for j=1:b
% counter1=counter1+1;
% if ISO(j)>=ST(j)+0.05 || ISO(j)>=ST(j)-0.05||ISO(j)==ST(j)
% counter2=counter2+1;
% else
% counter3=counter3+1;
% end
% 
% end
% end
% 
% % clear ISO;
% % clear ST;
% 
% % if counter3/counter1>=0.9
% %     warn = 'WARNING: ST elevation';
% %     fprintf(1,'\nK> WARNING: ST elevation \n');
% % else
% % %     warn = 'No MI';
% %    fprintf(1,'\nK> No MI \n');
% % end
% 
% 
% 
% % if counter3/counter1>=0.95 ||mean_ST_dvalue >1 ||mean_path_c<=-0.25||
% %     warn = 'Symptom of MI exists, need to check BP and Oxyzen saturation ';
% %     set(handles.text10, 'String', warn);
% % else
% % %     warn = 'No MI';
% %     set(handles.text10, 'String', 'No symptom of MI');
% % end
% 
% % if counter3/counter1>=0.9 ||mean_ST_dvalue <-.03 
% %        e=4;
% % else
% %     e=0;
% % end
% % if mean_path_c<=-0.25 ||H_T_peak >.5
% %     f=1;
% %     else
% %     f=0;
% % end
% % 
% % if  mean_path_QT >.4
% %     g=2;
% %     else
% %     g=0;
% % end 
% %     if (e+f+g)>=3 && ((c_bp >= 120 & c_bp < 160)||(c_bp <105 & c_bp >=90))  && ((c_ox <93 & c_ox >=88)|| H_R>80)
% %         
% %         warn1 = 'Medium case of MI ';
% %     fprintf(1,'\nK> Medium case of MI \n');
% %     elseif (e+f+g)>=3 && (((c_bp >= 120 & c_bp < 160)||(c_bp <105 & c_bp >=90)) || H_R>80 || (c_ox <93 & c_ox >=88))
% %         
% %         warn2 = 'Initial case of MI ';
% %     fprintf(1,'\nK> Initial case of MI  \n');   
% %     
% %     
% % elseif (e+f+g)>3 && (c_bp>160 || c_bp<90)&& c_ox<88
% %     warn3 = 'Severe case of MI ';
% %     fprintf(1,'\nK> Severe case of MI  \n');  
% % elseif ((e+f+g)>0 & (e+f+g)<3) && ( c_ox>93 || (c_bp >= 105 & c_bp <= 120))
% % % elseif file_bp >90
% % 
% %    
% %     fprintf(1,'\nK> Potential MI scenario , need monitoring  \n'); 
% % elseif ((e+f+g)>0 & (e+f+g)<3) && ((c_bp >= 120||c_bp<105 ) || H_R>80 || c_ox <93 )
% % % elseif file_bp >90
% % 
% %     
% %      fprintf(1,'\nK> Potential MI scenario , need monitoring \n'); 
% %     
% %     
% %     
% % elseif (e+f+g)>0 && (f+g)<3 && c_ox>93 && (c_bp >= 105 & c_bp <= 120)
% %   
% %      fprintf(1,'\nK> need ECG monitoring for consistancy \n'); 
% % elseif H_R>80 && (e+f+g)<=1 
% %     
% %      fprintf(1,'\nK> Arrhythmia detected , need further diagnostic for arrhythmia type \n'); 
% %     
% %     else
% % 
% %     fprintf(1,'\nK> No MI symptoms\n'); 
% %     end
% 
% 
% % figure 
% %  plot(t,QW,'b');hold on;
% % %    plot (P_peak_t,P_peak, '*m');hold on;
% %      plot (P_t,P_amp,'*m');
% %    plot(S_t,S_amp,'+r'), grid on; hold on;
% %  plot(R_t,R_amp,'+k'); hold on;
% %  xlabel('Time');ylabel('Amplitude');
% %     plot (Q_t, Q_amp, 'o'); hold on;
% %    plot (TP_t,TP_amp, '+m');hold on;
% %   plot (T_peak_t,T_peak, 'o');hold on;
% %   plot (J_t,J_amp,'+r');hold on;
% %     plot (K_t,K_amp,'*r');hold on;
% % 
% % title('ECG signal with PQRSTJK point')
% % ylabel('ECG+S+R+Q+P+J')
% % hold off;
% % 
% %  figure 
% % 
% % plot(t(1:R_index(2)+80),QW(1:R_index(2)+80),'b');hold on;
% %   
% %    plot(S_t(1:2),S_amp(1:2),'+r'), grid on; hold on;
% %  plot(R_t(1:2),R_amp(1:2),'+k'); hold on;
% %  xlabel('Time');ylabel('Amplitude');
% %     plot (Q_t(1:2), Q_amp(1:2), 'o'); hold on;
% % 
% % figure 
% % 
% % %    plot (P_peak_t(1:2),P_peak(1:2), '*m');hold on;
% %      plot (P_t(1:2),P_amp(1:2),'*m');
% % 
% %    
% %    plot (TP_t(1:2),TP_amp(1:2), '+m');hold on;
% %   plot (T_peak_t(1:2),T_peak(1:2), 'o');hold on;
% %   plot (J_t(1:2),J_amp(1:2),'+r');hold on;
% %     plot (K_t(1:2),K_amp(1:2),'*r');hold on;
% % 
% % title('ECG signal with PQRSTJK point')
% % ylabel('ECG+S+R+Q+P+J')
% % hold off;
% % figure 
% %  plot(t,QW,'b');hold on;
% %   
% %    plot(S_t,S_amp,'+r'), grid on; hold on;
% %  plot(R_t,R_amp,'+k'); hold on;
% %  xlabel('Time');ylabel('Amplitude');
% %     plot (Q_t, Q_amp, 'o'); hold on;
% %   
% % 
% % title('QRS complex ')
% % ylabel('Amplitude')
% % xlabel('time');
% % hold off;
% % 
% % 
% % figure 
% %  plot(t(1:R_index(2)+80),QW(1:R_index(2)+80),'b');hold on;
% %   
% %    plot(S_t(1:2),S_amp(1:2),'+r'), grid on; hold on;
% %  plot(R_t(1:2),R_amp(1:2),'+k'); hold on;
% %  xlabel('Time');ylabel('Amplitude');
% %     plot (Q_t(1:2), Q_amp(1:2), 'o'); hold on;
% %   
% % 
% % title('QRS complex ')
% % ylabel('Amplitude')
% % xlabel('time');
% % hold off;
% 
% 
% 
% end
% % H_R_PPG = round(60/(mean (HR_PPG)))  
% % H_R = round(60/(mean (HR)))
% % DP=min(D_BP)
% % SP=max(D_BP)
% 
