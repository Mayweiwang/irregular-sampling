
       clear; close all;
       colors={'-ro','-b+','-m^','-k*','-gx','-cd'};

%        figure(5)
       N=64;m=2;
       tao=0.39;fr=13.385;
       Fs = 1/tao;   
       A1=1; A2=1;
       vr=0.2;
       fc=2*vr*fr/300;
       fb=0.102*sqrt(fr);
       f1=-fb+fc; f2=fb+fc;
       sfs=43;sfi=24;

%% plot out C(?, ?)
       n=-Fs/2:Fs/(N-1):Fs/2;
       for i=1:5
           sigma=.25*2*(i-1)*tao;
           omega=n*2*pi;
           c=exp(-sigma^2/2*omega.^2);
           plot(0:1/(N/2-1):1, c(N/2+1:N),colors{i});hold on;
       end
        legend('\sigma = T/4','\sigma = T/2','\sigma = T'), grid on;
        title(' \fontsize{15} C(\omega, \sigma)');
        set(gcf,'Color',[1 1 1]);
        close all;
   
%%       
       y=zeros(4,N);sigmae=zeros(1,300);
       for k=1:300
           
                   ind=0:N-1;
                   t = ind*tao;
                   delt=normrnd(0,sigma,1,length(t));
                   x1=A1*exp(sqrt(-1)*(2*pi*t*f1))+A2*exp(sqrt(-1)*(2*pi*t*f2)); 
                   x2 = A1*exp(sqrt(-1)*(2*pi*(t+delt)*f1))+A2*exp(sqrt(-1)*(2*pi*(t+delt)*f2));  
        
                   h = spectrum.periodogram('rectangular');
                   hopts1 = psdopts(h,x1);  % Default options based on the signal x
                   set(hopts1,'Fs',Fs,'SpectrumType','twosided','CenterDC',true,'NFFT',N);
                   subplot(221), psd(h,x1,hopts1);
                   h1=psd(h,x1,hopts1); Pxx1 = h1.Data; 
                   y(1,:)=10*log10(Pxx1);

                   hopts2 = psdopts(h,x2);  % Default options based on the signal x
                   set(hopts2,'Fs',Fs,'SpectrumType','twosided','CenterDC',true,'NFFT',N);
                   subplot(223), psd(h,x2,hopts2)
                   h2=psd(h,x2,hopts2); Pxx2 = h2.Data; y(2,:)=10*log10(Pxx2);
                   set(gcf,'color',[1 1 1]);

                   ind = 0:N*8-1;
                   t = ind*tao;
                   delt=normrnd(0,sigma,1,length(t));
                   x3 = A1*exp(sqrt(-1)*(2*pi*(t+delt)*f1))+A2*exp(sqrt(-1)*(2*pi*(t+delt)*f2));  
                   h = spectrum.welch('rectangular',N,0);
                   hopts3 = psdopts(h,x3);  % Default options based on the signal x
                   set(hopts3,'Fs',Fs,'SpectrumType','twosided','CenterDC',true,'NFFT',N);
                   subplot(222), psd(h,x3,hopts3);
                   h3=psd(h,x3,hopts3); Pxx3 = h3.Data;    y(3,:)=10*log10(Pxx3);% useful
                   W=h3.Frequencies;                          %     plot(10*log10(Pxx3))
%                    sfs=167; sfi=93;  
                   sfs=42; sfi=23;                   
%                    x2ave=mean(x2fft,1);                             % amplitude average
                   sigma2=log(Pxx3(sfs)/Pxx3(sfi))/(4*pi^2*(W(sfi)^2-W(sfs)^2));
                   if sigma2 <0
                                 sigma2=-sigma2;
                   end
                   sigmae(k)=sqrt(sigma2); 
                   omega=W*2*pi;
                   c=exp(-sigmae(k)^2/2*omega.^2);
                   Pxx4=Pxx3./c.^2;
                   Pxx3(sfs)=Pxx4(sfs); Pxx3(sfi)=Pxx4(sfi); y(4,:)=10*log10(Pxx3);
                   subplot(224),hpsd=dspdata.psd(Pxx3,W,'Fs',Fs);plot(hpsd);
                   
                   
                   
   

         
                    ratio=y(:,sfs)-y(:,sfi);
                    tans4=10.^(ratio/2);
                    tans=tans4.^.25;
                    wdir(k,:)=atan(tans)*2/pi*180;  
          
       end

%        [f,x] = ksdensity(wdir(:,3));plot(x,f);
%        [f,xf] = ecdf(wdir(:,3));
%        ind1=find(xf<(wdir(1,1)+2));
%        ind2=find(xf(ind1)>(wdir(1,1)-2));
%        pdf2deg= f(ind2(end))-f(ind2(1));
      
       



%   for i=1:8
%                       delt=normrnd(0,sigma,1,N);
%                       ind2=ind(1:N)+(i-1)*N;
%                       x2 = A1*exp(sqrt(-1)*(2*pi*(t(ind2)+delt)*f1))+A2*exp(sqrt(-1)*(2*pi*(t(ind2)+delt)*f2));  
%                       x2fft(i,:)=fftshift(abs(fft(x2)));
%                       
%                    end



%      x2logave=mean(20*log10(x2fft),1);          % log average
%                    x2ave_c=x2ave./c;                                 % calibrated
%                     y(1,:)=20*log10(xo);                             
%                     y(2,:)=20*log10(xbad);
%                     y(3,:)=x2logave;
%                     y(4,:)=x2logave-20*log10(c);
%                     y(5,:)=20*log10(x2ave_c);
%                     figure(1)
%                     plot(y(1,:));
%                     plot(n,y(1,:),colors{1}); hold on;
%                     plot(n,y(2,:),colors{2}); 
%                     plot(n,y(3,:),colors{3}); 
%                     legend('regular sampled sea echo','sea echo with sampling jitters','emsemble log average'), grid on;
%                     title(' \fontsize{15} First-order sea echo power spectrum');
%                     set(gcf,'Color',[1 1 1]);
%                     figure(2)
%                     plot(n,y(1,:),colors{1}); hold on;
%                     plot(n,y(4,:),colors{4}); 
%                     plot(n,y(5,:),colors{5}); 
%                     legend('regular sampled sea echo','ensemble log average with calibration','emsemble amplitude average with calibration'), grid on;
%                     title(' \fontsize{15} First-order sea echo power spectrum');
%                     set(gcf,'Color',[1 1 1]);
                
%                     plot(n,y(2,:),'y'); pause; hold on; plot(n,y(3,:),'m');
%                     pause; plot(n,y(1,:),'k'); pause;
%                     plot(n,y(4,:),'r'), hold on; pause;plot(n,y(5,:),'b'); pause;  
       
       