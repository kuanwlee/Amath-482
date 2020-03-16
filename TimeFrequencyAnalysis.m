clear all 
close all
clc

load handel
v = y';
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

L = length(v)/Fs;
n = length(v);
t2 = linspace(0,L,n+1); t = t2(1:n);
k=(2*pi/(2*L))*[0:(n/2) (-n/2):-1]; ks=fftshift(k);

v_ft = fft(v);

%% Different width of Gabor window
tslide = 0:0.1:L;
a_vec = [1 2 5 0.1];
for jj = 1:length(a_vec)
    sgt_spec = zeros(length(tslide),n);
    a = a_vec(jj);
    for i = 1:length(tslide)
        g = exp(-a*(t-tslide(i)).^2);
        sg = v.*g;
        sg_ft = fft(sg);
        sgt_spec(i,:) = abs(fftshift(sg_ft));
    end    
    subplot(2,2,jj)
    pcolor(tslide,ks,sgt_spec.')
    shading interp
    title(['a = ',num2str(a)],'Fontsize',16)
    set(gca,'Ylim',[-3000 3000],'Fontsize',16)
    ylabel('\omega')
    xlabel('Time(s)')
    colormap(hot)
end

%% Different step size of Gabor window

stpsize = [0.1 1 0.8 0.01];
a = 1;
for jj = 1:length(stpsize)
    tslide = 0:stpsize(jj):L;
    sgt_spec = zeros(length(tslide),n);
    
    for i = 1:length(tslide)
        g = exp(-a*(t-tslide(i)).^2);
        sg = v.*g;
        sg_ft = fft(sg);
        sgt_spec(i,:) = abs(fftshift(sg_ft));
    end    
    subplot(2,2,jj)
    pcolor(tslide,ks,sgt_spec.')
    shading interp
    title(['stepsize = ',num2str(stpsize(jj))],'Fontsize',12)
    set(gca,'Ylim',[-8000 8000],'Fontsize',12)
    ylabel('\omega')
    xlabel('Time(s)')
    colormap(hot)
end

%% Different shape of Gabor window

tslide = 0:0.1:L;
sgt_spec = zeros(length(tslide),n);
a = 1;    
for i = 1:length(tslide)
    g = (1-(t-tslide(i)).^2).*exp(-a*(t-tslide(i)).^2);
    sg = v.*g;
    sg_ft = fft(sg);
    sgt_spec(i,:) = abs(fftshift(sg_ft));
end    
subplot(1,2,1)
pcolor(tslide,ks,sgt_spec.')
shading interp
title('Mexican hat ','Fontsize',12)
set(gca,'Ylim',[-8000 8000],'Fontsize',12)
ylabel('\omega')
xlabel('Time(s)')
colormap(hot)

for i = 1:length(tslide)
    a = 0.5;
    width = t-tslide(i);
    g = zeros(1,length(t));
    g(0<=width&width<=a/2) = 1;
    g(a/2<=width&width<=a) = -1;
    sg = v.*g;
    sg_ft = fft(sg);
    sgt_spec(i,:) = abs(fftshift(sg_ft));
end    
subplot(1,2,2)
pcolor(tslide,ks,sgt_spec.')
shading interp
title('Step function ','Fontsize',12)
set(gca,'Ylim',[-8000 8000],'Fontsize',12)
ylabel('\omega')
xlabel('Time(s)')
colormap(hot)

%% Part 2 of HW2
clear all;close all; clc
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');

L = length(y')/Fs;
n = length(y');
t2 = linspace(0,L,n+1); t = t2(1:n);
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
% p8 = audioplayer(y,Fs); playblocking(p8);
%%
tslide = 0:0.1:L;
sgt_spec = zeros(length(tslide),n);
N = length(y');    
dF = Fs/N;
f = -1*(Fs/2):dF:Fs/2-dF;
tau = 1;
filter = zeros(1,length(y'));
    a = 40;
    for i = 1:length(tslide)
        g = exp(-a*(t-tslide(i)).^2);
%         filter(f-500<=0) = 1;
        sg = y'.*g;
%         sg = sg.*filter;
        sg_ft = fft(sg);
%         sg_ft = fftshift(sg_ft).*filter;
        sgt_spec(i,:) = abs(fftshift(sg_ft));
%         sgt_spec(i,:) = abs(sg_ft);
    end  
N = length(y');    
dF = Fs/N;
f = -1*(Fs/2):dF:Fs/2-dF;

figure(2)

pcolor(tslide,f,sgt_spec.')
shading interp
set(gca,'Ylim',[0 1000],'Fontsize',12)
title('Mary had a little lamb (Piano)')
ylabel('frequency (HZ)')
xlabel('Time(s)')
colormap(hot)

%%

clear all;close all; clc
[y,Fs] = audioread('music2.wav');
tr_piano=length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');
%%
L = length(y')/Fs;
n = length(y');
t2 = linspace(0,L,n+1); t = t2(1:n);
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
p8 = audioplayer(y,Fs); playblocking(p8);
tslide = 0:0.1:L;
sgt_spec = zeros(length(tslide),n);
N = length(y');    
dF = Fs/N;
f = -1*(Fs/2):dF:Fs/2-dF;
tau = 1;
filter = zeros(1,length(y'));
    a = 40;
    for i = 1:length(tslide)
        g = exp(-a*(t-tslide(i)).^2);
        filter(f-500<=0) = 1;
        sg = y'.*g;
%         sg = sg.*filter;
        sg_ft = fft(sg);
%         sg_ft = fftshift(sg_ft).*filter;
        sgt_spec(i,:) = abs(fftshift(sg_ft));
%         sgt_spec(i,:) = abs(sg_ft);
    end  
N = length(y');    
dF = Fs/N;
f = -1*(Fs/2):dF:Fs/2-dF;

figure(2)

pcolor(tslide,f,sgt_spec.')
shading interp
set(gca,'Ylim',[0 2000],'Fontsize',12)
title('Mary had a little lamb (recorder)')
ylabel('frequency (HZ)')
xlabel('Time(s)')
colormap(hot)




 
