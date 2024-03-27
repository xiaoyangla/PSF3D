clc,clear,close
path = "D:\besselData1\";

%% param design for fraunhofer propagation

% 1/dleta_x1 >= 2*Bx
% Bx=D/2/lambda/z
% L2 = lambda*z/delta_x1
% delta_x2 = lambda*z/L1
% N*dx2/lambda/z=1/dx1
% 1/dx2 = fc ; Fourier plane frequency
%%
pixelsize = 0.056;   
 
N1 = 42001;           % pupil plane sampling numbers
N = 2001;              % Fourier plane sampling numbers
lambda = 0.52;        % wavelength,unit um
NA = 0.734;           % objective NA
r0 = 5.9;             % aperture radius max of objective;um
f = 8e3;              %focal length;unit um
z = f;                % distance from Z0

L1 = lambda*z/pixelsize;%side length;unt mm,general L=2(2a) or 3(2a)
dx1 = L1/N1;   % pupil plane samping interval
%  x1 = linspace(-L1/2,L1/2,N);% sampling interval ,unit mm
% dx11 = lambda*z/N1/pixelsize;
x1 = single(-L1/2:dx1:L1/2-dx1); 
%dx1 = x1(2)-x1(1);  % sample interval
%x1 = single(-L1/2:dx1:L1/2-dx1);
%[X1,Y1] = ndgrid(x1,y1);

%SampleIntensity = SampleIntensity ./ max(max(SampleIntensity)); 

k = 2*pi/lambda;     % wave vector

% Bx = 2*pi/tao or 1/tao; tao is 脉冲宽度;
% deltaFreq =2pi/T,Tis 信号周期；离散信号周期无穷大，频谱线连续
% fs is samping frequency = 1/dx1, fueq_coord = f s/N*k, N is
% sampling numbers,k = (0,1,2,3...N-1)

a = 1e3;               % inner radius,um
delta = 0.1e3;         % split width,um
Na_eff = (a+delta/2)/f;
%fx = x/(lambda*z);  % frequency var to space var
%Bx = 10/L1;          % effctive bandwidth,10/(2a); dx<1/(2*Bx);
rhoMax_img = (a+delta)*2/(2*lambda*f);% fs = 1/dx1 >= 2*rhoMax_img,unit /um,1/rhoMax_img is space coords
%c = X1.^2 + Y1.^2;
%u1 = zeros(size(X1));
contant = (N-1)/2-300;
center =  (N1+1)/2;
plotScale = center-contant:center+contant;

u1 = single(singleRingMask_XYZ(x1,x1,a,delta));
%u1 = circshift(u1,[1,1]);     % shift
%u1(b^2<c&&c<a^2)=1;           % construct aperture 
Isrc = abs(u1.^2);             % src irradiance
figure(1)
subplot(1,2,1,'position', [0.05 0.15 0.5 0.5]);
%imagesc(x1,y1,Isrc); colorbar
imagesc(x1(plotScale),x1(plotScale),Isrc(plotScale,plotScale)); colormap('gray');
axis square;
axis xy;
xlabel('x （um）'); ylabel('y （um）');
title('Aperture Plane Irradiance');
hold on
hline = refline([0 0]);
hline.Color = 'r';
hline.LineStyle = '- -';
hold off
subplot(1,2,2,'position', [0.55 0.15 0.5 0.5]);
plot(x1(plotScale),Isrc(center,plotScale),'r');
xlabel('x （um）'); ylabel('y （um）');
axis square;
axis xy;
title('X Profile');
saveas(1,strcat(path,'AperturePlaneIrradiance.jpeg'))
clear Isrc

% dx2 = pixelsize;
%fraunhofer diffraction field; fourier transform of aperture
%[uh,L2] = propFrahDiff(u1,L1,lambda,z,z);
%[uh] = propFrahDiff(u1,dx2,dx1,lambda,z,z);
[uh,L2,x2] = propFrahDiff(u1,dx1,L1,lambda,z,z);
%dx2 = L2/N; %
% L2=lambda*z/dx1; %obs sidelength
dx2=lambda*z/L1; %obs sample interval
clear u1

% L2 = dx2*N;
% dy2 = dx2;
% x2=-L2/2:dx2:L2/2-dx2; %obs ords

uh_amp = abs(uh);
Ih = uh_amp.^2;
Ih_norm = NormalizeArr(Ih,0,1);
%U1 = ifftshift(fft2(fftshift(u1)))/lambda/f*dx1*dy1; %shift,2D fft,dxdy scaling；FFT：{?...dxdy} ——>{∑∑...delxdely},but fft2 do not *dxdy
%U11 = abs(fftshift(U1));   %center,get amptitud

%Iobs_nor = NormalizeArr(Iobs,0,1);
figure(2)
subplot(1,2,1,'position', [0.05 0.15 0.5 0.5]);
imagesc(x2(plotScale),x2(plotScale),Ih_norm(plotScale,plotScale));colorbar 
%line(fx(N/2-49:N/2+50)*lambda*z,N/2+1,'Color','r','LineWidth',5);
%colormap('gray');
xlabel('x （um）'); ylabel('y （um）');
axis square;
axis xy;
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off
title('Ariy Disc');
subplot(1,2,2,'position', [0.55 0.15 0.5 0.5]);
plot(x2(plotScale),Ih_norm(center,plotScale),'r');
xlabel('x （um）'); ylabel('Magnitude');
axis square;
axis xy;
fwhm = calFWHM(x2(plotScale),Ih_norm(center,plotScale));
fwhm = strcat('FWHM = ',num2str(fwhm),'um');
title(fwhm)
%[~,~,PSFz_fwhm,~] =findpeaks(OverallCrossSection,'NPeaks',1,'SortStr','descend'); FWHM
saveas(2,strcat(path,'AriyDisc.jpeg'))

figure(3)
subplot(1,2,1,'position', [0.05 0.15 0.5 0.5]);
%s = 50;% display width
imagesc(x2(plotScale),x2(plotScale),uh_amp(plotScale,plotScale));colorbar 
%line(fx(N/2-49:N/2+50)*lambda*z,N/2+1,'Color','r','LineWidth',5);
%colormap('gray');
xlabel('x （um）'); ylabel('y （um）');
axis square;
axis xy;
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off
title('PSF');
subplot(1,2,2,'position', [0.55 0.15 0.5 0.5]);
plot(x2(plotScale),uh_amp(center,plotScale),'r');
xlabel('x （um）'); ylabel('Magnitude');
axis square;
axis xy;
saveas(3,strcat(path,'PSF.jpeg'))

clear uh_amp

% transform airy   function back
%% fftshift shift right N/2+1 to index 1
%  ifftshift shift left
Ih = Ih(plotScale,plotScale);

N2 = 105;
center1 = (length(plotScale)+1)/2;
plotScale1 = center1-N2:center1+N2;
if mod(N1,2)==0
    Ihf = ifftshift(fft2(fftshift(Ih)))*(dx2*dx2);
else
    Ihf = fftshift(fft2(ifftshift(Ih)))*(dx2*dx2);
end

Ihf_amp = abs(Ihf);
Ihf_amp_nor = NormalizeArr(Ihf_amp,0,1);

%freq coords
fx = -1/(2*dx2):1/dx2/length(plotScale):1/(2*dx2)-1/dx2/length(plotScale);%freq coords,fnx = 1/(2*dx),frequency sample intervals:1/L2,
%fx = fx(center-N2:N2+center);
figure(4)
subplot(1,2,1,'position', [0.05 0.15 0.5 0.5]);
imagesc(fx(plotScale1),fx(plotScale1),Ihf_amp_nor(plotScale1,plotScale1));  colorbar;% imaginary part is tiny
%imagesc(fx,fy,Ihf_amp_nor);  colorbar;% imaginary part is tiny
axis square;
axis xy;
xlabel('fx （cyc/um）'); ylabel('fy（cyc/um）');
title('OTF');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

subplot(1,2,2,'position', [0.55 0.15 0.5 0.5]);
plot(fx(plotScale1),Ihf_amp_nor(center1,plotScale1),'r');
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('fx （cyc/um）'); ylabel('magnitude');
saveas(4,strcat(path,'OTF.jpeg'))

clear Ih Ihf_amp_nor Ihf Ihf_amp Ih_norm x1
% compare fft result to geometrical overlap of two disks.
% for a slice through the center, abs(x) equals the radius
% censlice = zz((N+1)/2,:);
% ra = abs(x1/a);
% olap = (pi*a.^2)*(2/pi)*(acos(ra/2)-(ra/2).*sqrt(1-(ra/2).^2));
% figure(4)
% plot(x1,olap,x1,censlice,'.')

dx1 = 0.06e3;% unit um
x2 = x2(plotScale);
uh = uh(plotScale,plotScale);
deltaZ = 0.056*4;  %unit um
sampling = length(plotScale);
zd = (-(sampling-1)/2:(sampling-1)/2)*deltaZ;
center2 = (sampling+1)/2;
%%
PSF_3D = single(zeros(sampling,sampling,length(zd)));
for i=1:length(zd)
    PSF_3D(:,:,i) = single(abs(propFresDiff(uh,dx2,lambda,zd(i))));
end
close all
clear uh
PSF_3D = PSF_3D.^2;
name = strcat("illumination_PSF_OTF_3D_",'Na_eff_',strrep(num2str(Na_eff),'.','p'),'.mat');
% % %sprintf('%s/%s/decon_simulation/', outputdir, subfolder);
filename = strcat(path,name);
PSF_3D = NormalizeArr(PSF_3D,0,1);
% OTF_3D = abs(ifftshift(fftshift(fftn(PSF_3D))));
OTF_3D = abs(fftshift(fftn(ifftshift(PSF_3D))));
OTF_3D = NormalizeArr(OTF_3D,0,1);
save(filename,'PSF_3D','OTF_3D','-v7.3')
fprintf('%s has been saved',name)
% formatSpec = 'X is %4.2f meters or %8.3f mm\n';
% fprintf(formatSpec,A1,A2)

% %[minA,maxA] = bounds(PSF_3D,1);
YZ_PSF = squeeze(PSF_3D(center1,:,:));
XY_PSF = squeeze(PSF_3D(:,:,center2));
save(strcat(path,'PSF_2D_ill.mat'),'YZ_PSF','XY_PSF','-v7.3')
fprintf('%s has been saved','PSF_2D_ill.mat')
% %[minA,maxA] = bounds(maxA,2); 
% %maxA = squeeze(maxA);
%%
load(filename);
load(strcat(path,'PSF_2D_ill.mat'))
ill_YZ_PSF = YZ_PSF;
ill_XY_PSF = XY_PSF;
figure(5)
subplot(1,2,1,'position', [0.05 0.15 0.5 0.5]);

%s = N/20;
imagesc(zd,x2,ill_YZ_PSF);  colorbar;% imaginary part is tiny
%imagesc(fx,fy,Ihf_amp_nor);  colorbar;% imaginary part is tiny
axis square;
axis xy;
xlabel('z （um）'); ylabel('x（um）');
title('YZ-PSF-ill');
% hold on
% hline = refline([0 0]);%add reference line
% hline.Color = 'r';
% hline.LineStyle = '- -';
% hold off

subplot(1,2,2,'position', [0.55 0.15 0.5 0.5]);
%ss =100;
plot(zd,ill_YZ_PSF(center1,:),'r'); 
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('z （um）'); ylabel('magnitude');

fwhm = calFWHM(zd,ill_YZ_PSF(center1,:));
fwhm = strcat('FWHM = ',num2str(fwhm),'um');
title(fwhm)
saveas(5,'D:\besselData1\YZ_PSF.jpeg')

load('D:\besselData\pSF_det_px_56nm\xz_Det_PSF_OTF_norm_520_NA0p74_px_56nm_RichardsWolf.mat');
xz_Det_PSF = permute(squeeze(xz_PSF_RW_520nm_NA0p734),[2,1]);
[M,N] = size(xz_Det_PSF);
deltaSize = (length(plotScale)-M)/2;
center3 = (M+1)/2;
plotScale3 = center3-50:center3+50;

x3 = x2(deltaSize+1:end-deltaSize);
fx2 = fx(deltaSize+1:end-deltaSize);

overall_PSF = ill_XY_PSF(deltaSize+1:end-deltaSize,deltaSize+1:end-deltaSize).*xz_Det_PSF;
overall_PSF_norm = NormalizeArr(overall_PSF,0,1);
overall_OTF = abs(fftshift(fftn(ifftshift(overall_PSF_norm))));
overall_OTF_norm= NormalizeArr(overall_OTF,0,1);

figure(6)
subplot(1,3,1,'position', [0.05 0.15 0.3 0.5]);
imagesc(x3(plotScale3),x3(plotScale3),overall_PSF_norm(plotScale3,plotScale3));  colorbar;
axis square;
axis xy;
xlabel('x（um）'); ylabel('z（um）');
title('Overall PSF');

subplot(1,3,2,'position', [0.35 0.15 0.3 0.5]);
plot(x3(plotScale3),overall_PSF_norm(center3,plotScale3),'r'); 
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('z （um）'); ylabel('magnitude');
fwhm = calFWHM(x3(plotScale3),overall_PSF_norm(center3,plotScale3));
fwhm = strcat('FWHM = ',num2str(fwhm),'um');
title(fwhm)

subplot(1,3,3,'position', [0.65 0.15 0.3 0.5]);
plot(x3(plotScale3),overall_PSF_norm(plotScale3,center3),'r'); 
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('x （um）'); ylabel('magnitude');

fwhm = calFWHM(x3(plotScale3),overall_PSF_norm(plotScale3,center3));
fwhm = strcat('FWHM = ',num2str(fwhm),'um');
title(fwhm)

saveas(6,strcat(path,'overall_PSF.jpeg'))


figure(7)
subplot(1,3,1,'position', [0.05 0.15 0.3 0.5]);
imagesc(fx2(plotScale3),fx2(plotScale3),overall_OTF_norm);  colorbar;
axis square;
axis xy;
xlabel('fx （cyc/um）'); ylabel('fz（cyc/um）');
title('Overall OTF');

subplot(1,3,2,'position', [0.35 0.15 0.3 0.5]);
plot(fx2(301:701),overall_OTF_norm(center3,301:701),'r'); 
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('fx （cyc/um）'); ylabel('magnitude');

subplot(1,3,3,'position', [0.65 0.15 0.3 0.5]);
plot(fx2(plotScale3),overall_OTF_norm(plotScale3,center3),'r'); 
%plot(fx,Ihf_amp_nor); 
axis square;
axis xy;
xlabel('fz （cyc/um）'); ylabel('magnitude');
saveas(7,strcat(path,'overall_OTF_norm.jpeg'))



% DitheredPSF_ill = repmat(YZ_PSF(:,1001),1,1001);
% figure(6)
% imagesc(DitheredPSF_ill);  colorbar;% imaginary part is tiny
% %imagesc(fx,fy,Ihf_amp_nor);  colorbar;% imaginary part is tiny
% axis square;
% axis xy;
% xlabel('X'); ylabel('Z');
% title('Dithered bessel');



%close all
%% sanjiaohanshu
trirec = zeros(size(X1));
b = 0.002;
L3 = 5*b;
dx3 = L3/N;
dy3 = dx3;
x3 = -L3/2:dx3:L3/2-dx3;
y3 = x3;
[X3,Y3] = ndgrid(x3,y3);
trirec(abs(X3)+abs(Y3)<b)=1;% construct aperture     
s = 800;

figure(5)
subplot(1,3,1);
imagesc(x3(N/2-s:N/2+s),y3(N/2-s:N/2+s),trirec(N/2-s:N/2+s,N/2-s:N/2+s)); colormap('gray');
axis square;
axis xy;

s = 100;
subplot(1,3,2);
imagesc(x2(N/2-s:N/2+s),x2(N/2-s:N/2+s),Ih(N/2-s:N/2+s,N/2-s:N/2+s)); colormap('gray');
axis square;
axis xy;
%% conv2(a,b)
z2 = conv2(trirec,Ih,'same');% convolution    
subplot(1,3,3);
s = 800;
imagesc(x3(N/2-s:N/2+s),y3(N/2-s:N/2+s),z2(N/2-s:N/2+s,N/2-s:N/2+s)); colormap('gray');
axis square;
axis xy;
[z2f,fxx,fyy] = FFT2D(z2,dx3,dy3,L3);
Iz2f_amp = abs(z2f);
Iz2f_amp_nor = NormalizeArr(Iz2f_amp,0,1);

[trirecF,fxxx,fyyy] = FFT2D(trirec,dx3,dy3,L3);
trirecF_amp = abs(trirecF);
trirecF_amp_nor = NormalizeArr(trirecF_amp,0,1);

z4 = trirecF.*Ihf;
Iz4 = abs(z4);
Iz4_nor = NormalizeArr(Iz4,0,1);

z3 = deconvwnr(z2,Ih,0.1);  
%U3 = FFT(I2),trirec
% c = 1;
% Z3 = conj(U2).*Gxy./(U2.*conj(U2)+1);
% Z33 = abs(Z3);
% Z33_nor = NormalizeArr(Z33,0,1);
[Iz3f,fxx1,fyy1] = FFT2D(z3,dx3,dy3,L3);
% Iz3f = abs(Gxy1);
% uz3 = ifftshift(ifft2(Z3));
Iz3f_amp = abs(Iz3f);
Iz3f_amp_nor = NormalizeArr(Iz3f_amp,0,1);
% Iz3f_nor = NormalizeArr(Iz3,0,1);

figure(6)
subplot(2,3,1);
imagesc(x3(N/2-s:N/2+s),y3(N/2-s:N/2+s),z2(N/2-s:N/2+s,N/2-s:N/2+s)); colormap('gray');
axis square;
axis xy;

subplot(2,3,4);
imagesc(x3(N/2-s:N/2+s),y3(N/2-s:N/2+s),z3(N/2-s:N/2+s,N/2-s:N/2+s)); colormap('gray');
axis square;
%axis xy;

s = 80 ;
subplot(2,3,2);
imagesc(fxxx(N/2-s:N/2+s),fyy(N/2-s:N/2+s),trirecF_amp(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
axis square;
%axis xy;
%xlabel('fx （cyc/mm）'); ylabel('fy（cyc/mm）');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

subplot(2,3,5);
%imagesc(fxx1(N/2-s:N/2+s),fyy1(N/2-s:N/2+s),Z33(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
imagesc(fxx1(N/2-s:N/2+s),fyy1(N/2-s:N/2+s),Iz3f_amp(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
axis square;
%axis xy;
%xlabel('fx （cyc/mm）'); ylabel('fy（cyc/mm）');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

s = 20 ;
subplot(2,3,3);
imagesc(fxxx(N/2-s:N/2+s),fyyy(N/2-s:N/2+s),Iz4_nor(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
%axis xy;
%xlabel('fx （cyc/mm）'); ylabel('magnitude');


subplot(2,3,6);
plot(fx(N/2+1:N/2+30*s),Ihf_amp_nor(N/2+1,N/2+1:N/2+30*s));
hold on;
plot(fxxx(N/2+1:N/2+s),trirecF_amp_nor(N/2+1,N/2+1:N/2+s));
hold on;
plot(fxxx(N/2+1:N/2+s),Iz4_nor(N/2+1,N/2+1:N/2+s));
hold on;
plot(fxx(N/2+1:N/2+s),Iz2f_amp_nor(N/2+1,N/2+1:N/2+s));
hold on;
plot(fxx1(N/2+1:N/2+s),Iz3f_amp_nor(N/2+1,N/2+1:N/2+s));
legend('hF','trirecF','trirecF.*Ihf','convolution','deconvolution')
axis square;
axis xy;
%xlabel('fx （cyc/mm）'); ylabel('magnitude');


%%
% proporgation, get Fz, angular spectra at d
%Fz = w22.*exp(1i*2*pi*d/lamda).*sqrt(1-(lamda.*).^2-(lamda.*).^2);
dx1 = 0.06;% unit mm
[u4] = propFresDiff(uh,L2,lambda,dx1);
I4 = abs(u4).^2;

figure(7)
subplot(3,3,1);
s=50;
imagesc(x2(N/2-s:N/2+s),x2(N/2-s:N/2+s),Ih(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('x （mm）'); ylabel('y（mm）');
title('Power');
subplot(3,3,2);
imagesc(x2(N/2-s:N/2+s),x2(N/2-s:N/2+s),I4(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('x （mm）'); ylabel('y（mm）');
title('Power');
I42 = I4+Ih;
subplot(3,3,3);
imagesc(x2(N/2-s:N/2+s),x2(N/2-s:N/2+s),I42(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('x （mm）'); ylabel('y（mm）');
title('Power');

subplot(3,3,4);
s = 600 ;
imagesc(fx(N/2-s:N/2+s),fy(N/2-s:N/2+s),Ihf_amp(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('fy（cyc/mm）');
title('FFT');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

U4 = fftshift(fft2(I4))*(dx2*dy2);
I4f = abs(U4).^2;
I4f_nor = NormalizeArr(I4f,0,1);

subplot(3,3,5);
s = 600 ;
imagesc(fx(N/2-s:N/2+s),fy(N/2-s:N/2+s),I4f(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('fy（cyc/mm）');
title('FFT');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

U42 = fftshift(fft2(I42))*(dx2*dy2);
I42f = abs(U42).^2;
I42f_nor = NormalizeArr(I42f,0,1);

subplot(3,3,6);
s = 700 ;
imagesc(fx(N/2-s:N/2+s),fy(N/2-s:N/2+s),I42f_nor(N/2-s:N/2+s,N/2-s:N/2+s));  colorbar;% imaginary part is tiny
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('fy（cyc/mm）');
title('FFT');
hold on
hline = refline([0 0]);%add reference line
hline.Color = 'r';
hline.LineStyle = '- -';
hold off

subplot(3,3,7);
plot(fx(N/2-s:N/2+s),Ihf_amp_nor(N/2,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('magnitude');

subplot(3,3,8);
plot(fx(N/2-s:N/2+s),I4f_nor(N/2,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('magnitude');

subplot(3,3,9);
plot(fx(N/2-s:N/2+s),I42f_nor(N/2,N/2-s:N/2+s));  colorbar;
axis square;
axis xy;
xlabel('fx （cyc/mm）'); ylabel('magnitude');

function[u2,L2,x2] = propFrahDiff(u1,dx1,L1,lambda,d1,z)
% propagation - Fraunhofer pattern
% assumes uniform sampling
% u1 - source plane field
% L1 - source plane side length
% lambda - wavelength
% z - propagation distance
% L2 - observation plane side length
% u2 - observation plane field
[M,N]=size(u1); %get input field array size
%dx1=L1/M; %source sample interval
k=2*pi/lambda; %wavenumber
%
L2=lambda*z/dx1; %obs sidelength
dx2=lambda*z/L1; %obs sample interval
L2=dx2*N;
x2=-L2/2:dx2:L2/2-dx2; %obs coords
%[X2,Y2]=meshgrid(x2,x2);
%
%pha = exp(1i*k/(2*z)*(1-d1/z)*(X2.^2+Y2.^2));
pha=1;
c=1/(1i*lambda*z)*pha;
if mod(M,2)==0
    u2=c.*ifftshift(fft2(fftshift(u1)))*dx1^2;
else
    u2=c.*fftshift(fft2(ifftshift(u1)))*dx1^2;
end


end

function[u2]=propFresDiff(u1,dx,lambda,z)
% propagation - transfer function approach
% assumes same x and y side lengths and
% uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field
[M,N]=size(u1); %get input field array size
%dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber
fx=-1/(2*dx):1/dx/N:1/(2*dx)-1/dx/N; %freq coords
[FX,FY]=meshgrid(fx,fx);
H=exp(1i*k*z)*exp(-1i*pi*lambda*z*(FX.^2+FY.^2)); %trans func
if mod(M,2)==0 
    H=fftshift(H); %shift trans func
else
    H=ifftshift(H); %shift trans func
end

if mod(M,2)==0
    U1=fft2(fftshift(u1)); %shift, fft src field
else
    U1=fft2(ifftshift(u1)); %shift, fft src field
end

U2=H.*U1; %multiply

if mod(M,2)==0
    u2=ifftshift(ifft2(U2)); %inv fft, center obs field
else
    u2=fftshift(ifft2(U2)); %inv fft, center obs field
end

end

function[u2]=propFresIR(u1,L,lambda,z)
% propagation - impulse response approach
% assumes same x and y side lengths and
% uniform sampling
% u1 - source plane field
% L - source and observation plane side length
% lambda - wavelength
% z - propagation distance
% u2 - observation plane field
[M,N]=size(u1); %get input field array size
dx=L/M; %sample interval
k=2*pi/lambda; %wavenumber
x=-L/2:dx:L/2-dx; %spatial coords
[X,Y]=meshgrid(x,x);
h=1/(1i*lambda*z)*exp(1i*k/(2*z)*(X.^2+Y.^2)); %impulse
H=ifftshift(fft2(fftshift(h)))*dx^2; %create trans func
U1=ifftshift(fft2(fftshift(u1))); %shift, fft src field
U2=H.*U1; %multiply
u2=ifftshift(ifft2(U2)); %inv fft, center obs field
end

function[Gxy,fx,fy] = FFT2D(gxy,dx,dy,L)

%g0 = fftshift(gxy);
G0 = fft2(gxy)*dx*dy;
Gxy = fftshift(G0);
fx = -1/(2*dx):1/L:1/(2*dx)-1/L;
fy = -1/(2*dy):1/L:1/(2*dy)-1/L;%freq coords

end

function[fwhm] = calFWHM(x,input)
input = squeeze(input);
% M = length(input);
%   if mod(M,2)==1
%      center1 = (M+1)/2;
%   else
%      center1 = M/2;
%   end
a = input;
[a_max,a_max_idx]= max(a);
idx1 = find(a(1:a_max_idx) <= a_max/2,1,'last');
idx2 = a_max_idx + find(a(a_max_idx:end) <= a_max/2,1,'first')-1;
fwhm = x(idx2)-x(idx1);
end