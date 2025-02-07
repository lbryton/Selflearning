function hh=example_nyq_sqrt(f_sym,f_smpl,alpha,m_dly,fig)

% function h=nyq_sqrt(f_sym,f_smpl,alpha,m_dly,fig);...
% f_sym  = symbol rate (nominally 1)
% f_smpl = samples/per symbol
% alpha  = Rolloff or excess BW (0< alpha < 1)
% m_dly  = number of symbols from start to filter peak (and center)
% fig    = 0 for no figure, 
%        = 1 for figure,.... Impulse Response and Frequency Response
% try y=nyq_sqrt(1,8,0.25,10,0);
% 
% written by fred harris, San Diego State University, 21-Sept 2013


M=f_smpl;
nn=-M*m_dly:M*m_dly;
n_ovr_M=nn/M;

z1=find(nn==0);                                  % denominator Zero 
z2=find(abs(1-(4*alpha*n_ovr_M).^2)<10^-10);     % denominator Zeros 
n2=z2-z1;

ZZ1=(4*alpha+pi*(1-alpha))/pi;              % 0/0 at z1
ZZ2=(-4*alpha*n2/M).*(pi*(1+alpha)/M).*sin(pi*(1+alpha).*n2/M)+ ...
     (4*alpha/M)*cos(pi*(1+alpha).*n2/M)+(pi*(1-alpha)/M)*cos(pi*(1-alpha).*n2/M);
 dd=(pi/M)*((1-(4*alpha.*n2/M).^2)-2*(4*alpha.*n2/M).^2);
ZZ2=ZZ2/dd;                                     % zeros at z2

tp=4*alpha*n_ovr_M.*cos(pi*(1+alpha)*n_ovr_M)+sin(pi*(1-alpha)*n_ovr_M);
bt=(1-(4*alpha*n_ovr_M).^2)*pi.*n_ovr_M;

   hh=(tp./bt);
   hh(z1)=ZZ1;
   hh(z2)=ZZ2;
   hh=hh/max(hh);
   
   if fig==1
    figure(1)
    subplot(2,1,1)
     plot(nn,hh,'.-')
     grid on
     axis([nn(1) -nn(1) -0.3 1.2])
     title('Impulse Response, SQRT Nyquist Filter')
     xlabel('Time Index')
     ylabel('Amplitude')
     text(m_dly/2,0.9,['alpha = ',num2str(alpha)])
    

    
     subplot(2,1,2)
     plot((-0.5:1/2048:0.5-1/2048)*f_smpl,fftshift(20*log10(abs(fft(hh/sum(hh),2048)))))
     hold on
     plot([-0.5 -0.5 0.5 0.5],[-80 0 0 -80],'r')
     hold off
     grid on
     axis([-f_smpl/2 f_smpl/2 -80 10])
     title('Frequency Response')
     xlabel('Normalized Frequency')
     ylabel('Log Mag (dB)')
   end