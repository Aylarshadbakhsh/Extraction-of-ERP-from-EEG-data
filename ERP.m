clc
clear
close all
%% load data
load('ex2data.mat');
%% sweeps numbers
T_1=length(indf);
T_2=length(indd);
%% time window (ms) -50-500 ms window
fs=250;
t=linspace(-50,500,137); %n=550ms*250/1000ms=137
%% part a 
y_e_1=average(eeg,indf);
y_e_2=average(eeg,indd);
figure()
subplot(2,1,1)
plot(t,y_e_1)
title('Regular stimuli')
xlabel('time')
ylabel('amplitude')
subplot(2,1,2)
plot(t,y_e_2)
title('Irregular stimuli')
xlabel('time')
ylabel('amplitude')
%% part b average of odd indexed elements and even indexed  elements
y_e_1even=average1even(eeg,indf);
y_e_1odd=average1odd(eeg,indf);
noise=(y_e_1odd-y_e_1even);
figure()
subplot(3,1,1 )
plot(t,y_e_1odd)
title('odd')
ylim([-4 3])
xlabel('time')
ylabel('amplitude')
subplot(3,1,2)
plot(t,y_e_1even)
title('even')
xlabel('time')
ylabel('amplitude')
subplot(3,1,3)
plot(t,noise)
title('noise')
xlabel('time')
ylabel('amplitude')
varnoise=var(noise);
varsig=var(y_e_1);
snr=varsig/varnoise;
fprintf(' SNR=%d.\n',snr);
 %% part b window(50 200)ms and  different sweeps
t2=linspace(15,200,46);
n1=10;
ERP2=average2(n1,indf,eeg);
ERP2odd=average2odd(n1,indf,eeg);
ERP2even=average2even(n1,indf,eeg);
figure()
plot(t2,ERP2)
title('n=10')
noise=ERP2odd-ERP2even;
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=10:%d.\n',snr);
n2=50;
ERP2=average2(n2,indf,eeg);
ERP2odd=average2odd(n2,indf,eeg);
ERP2even=average2even(n2,indf,eeg);
noise=ERP2odd-ERP2even;
figure()
plot(t2,ERP2)
title(' n=50 ')
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=50:%d.\n',snr);
n3=100;
ERP2=average2(n3,indf,eeg);
ERP2odd=average2odd(n3,indf,eeg);
ERP2even=average2even(n3,indf,eeg);
noise=ERP2odd-ERP2even;
figure()
plot(t2,ERP2)
title('n=100')
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=100:%d.\n',snr);
n4=200;
ERP2=average2(n4,indf,eeg);
ERP2odd=average2odd(n4,indf,eeg);
ERP2even=average2even(n4,indf,eeg);
noise=ERP2odd-ERP2even;
figure()
plot(t2,ERP2)
title('n=200')
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=200:%d.\n',snr);
n5=300;
ERP2=average2(n5,indf,eeg);
ERP2odd=average2odd(n5,indf,eeg);
ERP2even=average2even(n5,indf,eeg);
noise=ERP2odd-ERP2even;
figure()
plot(t2,ERP2)
title('n=300')
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=300:%d.\n',snr);
n6=400;
ERP2=average2(n6,indf,eeg);
ERP2odd=average2odd(n6,indf,eeg);
ERP2even=average2even(n6,indf,eeg);
figure()
plot(t2,ERP2)
title('n=400')
noise=ERP2odd-ERP2even;
varnoise=var(noise);
varsig=var(ERP2);
snr=varsig/varnoise;
fprintf(' SNR for n=400:%d.\n',snr);
n7=420;
y2=average2(n7,indf,eeg);
ERP2odd=average2odd(n7,indf,eeg);
ERP2even=average2even(n7,indf,eeg);
noise=ERP2odd-ERP2even;
figure()
plot(t2,y2)
title('n=420')
varnoise=var(noise);
varsig=var(y2);
snr=varsig/varnoise;
fprintf(' SNR for n=420:%d.\n',snr);

 %% part c FIR filter and IIR filter
%FIR
Fcp_low=1;  %lower cutoff frequency
Fcp_high=20;  %higher cutoff frequency
filterbandpass=fir1(7,[Fcp_low Fcp_high]/(250/2),'bandpass');
filter2=filtfilt(filterbandpass,1,y2);
figure()
plot(t2,filter2)
title('FIR')
xlabel('time')
ylabel('amplitude')
[h,w]=freqz(filter2,1,1000);
figure()
plot(w*(250/2),abs(h))
xlabel('frequency')
ylabel('magnitude')
title('frequency response of FIR filter')
%IIR
[z,p,k]=butter(7,[Fcp_low Fcp_high]/(250/2),'bandpass');
[sos,g]=zp2sos(z,p,k);
filter2=filtfilt(sos,g,y2);
figure()
plot(t2,filter2)
title('IIR')
xlabel('time')
ylabel('amplitude')
[h,w]=freqz(filter2,1,1000);
figure()
plot(w*250/2,abs(h))
title('frequency response IIR filter')
xlabel('frequency')
ylabel('magnitude')
 %% part c filtering eeg
Fcp_low=1;  %lower cutoff frequency
Fcp_high=20;  %higher cutoff frequency
[z,p,k]=butter(7,[Fcp_low Fcp_high]/(250/2),'bandpass');
[sos,g]=zp2sos(z,p,k);
filter2=filtfilt(sos,g,eeg);
n=T_1;
y1prime=average2(n,indf,filter2);
nprime=T_2;
y2prime=average2(nprime,indd,filter2);
subplot(2,1,1)
plot(t2,y1prime)
title(' filtered response to regular stimuli')
xlabel('time')
ylabel('amplitude')
subplot(2,1,2)
plot(t2,y2prime)
xlabel('time')
ylabel('amplitude')
title('filtered response to irregular stimuli')
%% function average 
function y_e_1=average(eeg,s)
    T_1=length(s);
    ERP1=zeros(T_1,137);
    for i=1:T_1
     ERP1(i,:)=eeg(s(i)-12:s(i)+124);
    end
    mean_erp1=mean(ERP1.'); 
    for g=1:T_1
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end
%% function b average of odd indexed elements and even elements
function y_e_1=average1odd(eeg,s)
    T_1=length(s);
    ERP1=zeros(T_1,137);
    for i=1:2:T_1
     ERP1(i,:)=eeg(s(i)-12:s(i)+124);
    end
    mean_erp1=mean(ERP1.'); 
    for g=1:2:T_1
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end
function y_e_1=average1even(eeg,s)
    T_1=length(s);
    ERP1=zeros(T_1,137);
    for i=2:2:T_1
     ERP1(i,:)=eeg(s(i)-12:s(i)+124);
    end
    mean_erp1=mean(ERP1.'); 
    for g=2:2:T_1
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end

%% function part b (15 200)ms
function y_e_1= average2(n,s,eeg)
    ERP1=zeros(n,46);
    for i=1:n
     ERP1(i,:)=eeg(s(i)+8:s(i)+53);
    end
    mean_erp1=mean(ERP1.'); 
    for g=1:n
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end
function y_e_1= average2odd(n,s,eeg)
    ERP1=zeros(n,46);
    for i=1:2:n
     ERP1(i,:)=eeg(s(i)+8:s(i)+53);
    end
    mean_erp1=mean(ERP1.'); 
    for g=1:2:n
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end
function y_e_1= average2even(n,s,eeg)
    ERP1=zeros(n,46);
    for i=2:2:n
     ERP1(i,:)=eeg(s(i)+8:s(i)+53);
    end
    mean_erp1=mean(ERP1.'); 
    for g=2:2:n
     ERP1(g,:) = ERP1(g,:) - mean_erp1(g); 
    end
   y_e_1=mean(ERP1);
end