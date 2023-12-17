clear all; close all; clc;

[x,fs]=audioread("recenica33.wav");

peff=rms(x);
p0=2*10^(-5);
L=20*log10(peff/p0); %vrednost za ceo signal
N=length(x);
tau=0.03;
Todb=ceil(fs*tau);
Ts=1/fs;
xdc=x-mean(x); %uklanjamo DC komponentu
xn = xdc./max(abs(xdc)); %normalizujemo signal da bude u intervalu [-1 1]
M=tau/fs;
d=N/fs; %duzina trajanja signala u vremenu
Nodb=ceil(d/tau);
tosa=(0:N-1)/fs;

figure(1), plot(tosa, xn); title('vremenski oblik signala'); xlabel('t[s]');

%proracun za nivo signala
k=1;
for br= 1 : Todb/2 : N-Todb/2 %preklapanje 50%
    x1 =xn(br:br+Todb/2);
    L1(k)=20*log10(rms(x1)/p0);
    for br1 =2:Todb/2
        z1=abs(sign(x1(br1))-sign(x1(br1-1))); %pokusaj, koristila sam drugi metod
    end
    %zcr(k)=(z1/(Todb) )*M;
    [rate(k),zcr(k)]=zerocrossrate(x1); 
    k=k+1;
end

t=((0:k-2)*Todb)/(2*fs);

figure(2), plot(t,L1); xlabel('t[s]'); ylabel('[dB]'); title('nivo signala u dB');
figure(3), plotyy( tosa,xn,t,L1/max(L1)); title('vremenski oblik signala i intenzitet'); xlabel('t[s]');

figure(4), plot(t,zcr); xlabel('t[s]'); title('ZCR');
figure(5), plotyy(tosa,xn,t,zcr); title('vremenski oblik signala i ZCR'); xlabel('t[s]');

%zvucni i bezvucni
L_prag= 76;
zcr_prag=35;
odluka=ones(1,281);
for br= 1 : 281
    if((L1(br)<=L_prag)&&(zcr(br)>=zcr_prag))
    odluka(br)=0;
    end
end

figure(6), plotyy(t,odluka, tosa,xn); title('vremenski oblik signala i podela na zvucne i bezvucne foneme'); xlabel('t[s]'); legend( '0 - bezvucni, 1- zvucni');

%trazenje osnovne ucestanosti
f0=zeros(1,281);
for br= 1 : 281
    x1 =xn(br:br+Todb/2);
    if(odluka(br)==1)
        r=xcorr(x1,x1);
        [pks,loc] = findpeaks(r, 'SortStr', 'descend', 'NPeaks', 2); %nalazi 2 maksimuma
        f0(br)=abs(loc(1)-loc(2))/2;
    end
end
figure(7), plot(t,f0); xlabel('t[s]'); 
figure(8), plotyy( tosa,xn,t,f0); title('vremenski oblik signala i osnovna frekvencija zvucnih fonema'); xlabel('t[s]');

%trazenje formantnih oblasti
x1s=xn(1:ceil(N/d));%izdvajamo 1s
X1s=abs(fft(x1s));
X1s1=X1s(1:1000);
X1s2=X1s(1000:2000);
X1s3=X1s(2000:3000);
X1s4=X1s(3000:4800);
%na prezentaciji pise da govor ima jednu formantnu frekvenciju na svakih
%1000Hz pa sam zato pretrazivala maksimume na opsegu od 100Hz i priblizno
%je tacno
[f1,F1]=max(X1s1);
[f2,F2]=max(X1s2);
[f3,F3]=max(X1s3);
[f4,F4]=max(X1s4);
%F1 je 504Hz F2 je 1774 Hz F3 je 2683 F4 je 3864 iz Praata
[fpks,floc] = findpeaks(X1s, 'SortStr', 'descend', 'NPeaks', 4);
figure(9), plot(smooth(X1s));xlim([0 4000]); title('spektar signala', 'peakovi su frekvencije formantnih oblasti');
F2=F2+1000;
F3=F3+2000;
F4=F4+3000;
figure(10), plot(X1s);
hold on
plot(F1, f1, 'o');
hold on
plot(F2, f2, 'o');
hold on
plot(F3, f3, 'o');
hold on
plot(F4, f4, 'o');
xlim([0 5000]);title('spektar signala', 'peakovi su frekvencije formantnih oblasti');
%spektogram
[S,F,T,P]=spectrogram(x1s,512,512/2,fs,fs);
figure(11),
surf(T,F,10*log10(P),'edgecolor','none') 
view(0,90);
xlabel('Vreme (s)'); ylabel('Frekvencija (Hz)'); ylim([0 4000]); title('prva 4 formanta na spektogramu');
