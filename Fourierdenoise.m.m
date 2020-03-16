clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
Un_ft = zeros(length(Undata(:,1)),n,n,n);
for i = 1:length(Undata(:,1))
    Un(:,:,:) = reshape(Undata(i,:),n,n,n); 
    Un_ft(i,:,:,:) = fftn(Un(:,:,:),[n,n,n]);
end    

% isosurface(X,Y,Z,abs(Un),0.4)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
Un_ftavg = mean(Un_ft,1);
Un_ftavg3d(:,:,:) = Un_ftavg(1,:,:,:);
Un_ftavg3d = abs(Un_ftavg3d);
Un_ftavg3d = Un_ftavg3d./max(Un_ftavg3d,[],'all');
fv_fft = isosurface(Kx,Ky,Kz,fftshift(Un_ftavg3d),0.9);
isosurface(Kx,Ky,Kz,fftshift(Un_ftavg3d),0.9)
xlabel('kx')
ylabel('ky')
zlabel('kz')

kx = fv_fft.vertices(:,1);
ky = fv_fft.vertices(:,2);
kz = fv_fft.vertices(:,3);
tau = 0.2;
x_filt = mean(kx);
y_filt = mean(ky);
z_filt = mean(kz);
filter = exp(-1*tau*(Kx-(x_filt)).^2).*exp(-1*tau*(Ky-y_filt).^2).*exp(-1*tau*(Kz-(z_filt)).^2);
Un_new = zeros(20,n,n,n);
for m = 1:20
    Un_ft1(:,:,:) = Un_ft(m,:,:,:);
    Un_ft1new = fftshift(Un_ft1).*filter;
    Un_new(m,:,:,:) = ifftn(Un_ft1new,[n,n,n]);

end    

figure(2)
for i = 1:20
    Un_s(:,:,:) = Un_new(i,:,:,:);
    isosurface(X,Y,Z,abs(Un_s(:,:,:)),0.4)
    axis([-20 20 -20 20 -20 20]), grid on, hold on
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pause(0.1)
end
fv = isosurface(X,Y,Z,abs(Un_s(:,:,:)),0.4);

x_marb = mean(fv.vertices(:,1));
y_marb = mean(fv.vertices(:,2));
z_marb = mean(fv.vertices(:,3));
