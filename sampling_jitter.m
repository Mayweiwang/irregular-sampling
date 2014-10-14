
       clear; close all;
       

%        figure(5)
       N=64;m=2;
       tao=0.39;
       Fs = 1/tao;   
       A1=2; A2=.5;
       fb=0.102*sqrt(13.4);
       f1=-fb; f2=fb;
       sfs=42;sfi=24;
       sigma=1.5*tao;
       sigma0=1.5*tao;
       n=-Fs/2:Fs/(N-1):Fs/2;
       
       omega=n*2*pi;
       cb=exp(-sigma^2/2*omega.^2);
       cn=exp(-sigma0^2/2*omega.^2);
%        [pr sfs]=min(abs(n-fb));       [pl sfi]=min(abs(n+fb));
       
       
           
       c=cn;
       c(sfs-m:sfs+m)=cb(sfs-m:sfs+m); c(sfi-m:sfi+m)=cb(sfi-m:sfi+m);
     
       a=max([c(sfs) c(sfi)]);b=min([c(sfs) c(sfi)]);
       if A1 > A2
           c(sfs)=a;
           c(sfi)=b;
       else
             c(sfs)=b;
             c(sfi)=a;  
       end
      
       
       diff=zeros(300,4);diff2=zeros(300,4);y=zeros(5,N);
       for k=1:300
           
       for j=1:2
            if j==1
                   ind=1:N;
                   t = ind*tao;
                   delt=normrnd(0,sigma,1,length(t));
                   x1=A1*exp(sqrt(-1)*(2*pi*t*f1))+A2*exp(sqrt(-1)*(2*pi*t*f2)); 
                   x2 = A1*exp(sqrt(-1)*(2*pi*(t+delt)*f1))+A2*exp(sqrt(-1)*(2*pi*(t+delt)*f2));  

                   xo=fftshift(abs(fft(x1))); 
                   xbad=fftshift(abs(fft(x2)));
%                 plot(n,20*log10(xo),'r'); hold on;pause;
%                 plot(n,20*log10(xbad),'k');   pause
           else
                   ind = 1:64*8;
                   t = ind*tao;
                   delt=normrnd(0,sigma,1,length(t));
                   xfft=zeros(8,64);
                   x2fft=zeros(8,64);

                   for i=1:8
                      ind2=ind(1:64)+(i-1)*64;
                      x2 = A1*exp(sqrt(-1)*(2*pi*(t(ind2)+delt(ind2))*f1))+A2*exp(sqrt(-1)*(2*pi*(t(ind2)+delt(ind2))*f2));  
                      x2fft(i,:)=fftshift(abs(fft(x2)));
%                    plot(n,x2fft);
%                    title('no shift absolute complex number spectrum: exp');
                   end
                                      
                   x2ave=mean(x2fft,1);
                   x2logave=mean(20*log10(x2fft),1);

                    x2ave_c=x2ave./c;
                    y(1,:)=20*log10(xo);
                    y(2,:)=20*log10(xbad);
                    y(3,:)=x2logave;
                    y(4,:)=x2logave-20*log10(c);
                    y(5,:)=20*log10(x2ave_c);
%                     plot(n,y(2,:),'y'); pause; hold on; plot(n,y(3,:),'m');
%                     pause; plot(n,y(1,:),'k'); pause;
%                     plot(n,y(4,:),'r'), hold on; pause;plot(n,y(5,:),'b'); pause;

                    diff(k,:)=[y(1,sfs)-y(2,sfs) y(1,sfs)-y(3,sfs) y(1,sfs)-y(4,sfs) y(1,sfs)-y(5,sfs)];
                    diff2(k,:)=[y(1,sfi)-y(2,sfi) y(1,sfi)-y(3,sfi) y(1,sfi)-y(4,sfi) y(1,sfi)-y(5,sfi)];
                    ratio=y(:,sfs)-y(:,sfi);
                    tans4=10.^(ratio/2);
                    tans=tans4.^.25;
                    wdir(k,:)=atan(tans)*2/pi*180;  
           end
       end
       end

%        [f,x] = ksdensity(wdir(:,3));plot(x,f);
%        [f,xf] = ecdf(wdir(:,3));
%        ind1=find(xf<(wdir(1,1)+2));
%        ind2=find(xf(ind1)>(wdir(1,1)-2));
%        pdf2deg= f(ind2(end))-f(ind2(1));
      
       
       
       
       