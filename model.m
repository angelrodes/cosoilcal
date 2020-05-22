function [ C ] = model( P,L,l,z,effectivedepth,age,erosion,zs)
 % Calculates the concentration at several depths (zs) considering the
 % surface production rate data (P,L,l), a variable-density profile
 % (z,effectivedepth) for an age and several erosion rates
 
 % resize paramanters
 zs=zs(:); % vertical
 erosion=erosion(:)'; % horizontal
 z=z(:)'; % horizontal
 effectivedepth=effectivedepth(:)'; % horizontal
 age=min(age,4.5430e+09); % the surface can't be older than the Earth
 l=l+log(2)/(100*4543e6); % add a negligible decay (100 times earth age) to avoid 1/0
 
 
 % start concentration
 C=zeros(length(zs),length(erosion));
 
 time=logspace(2,log10(age),100);

 for n=1:length(time)-1
     Z=erosion*time(n)+zs;
     Zprev=erosion*time(n+1)+zs;
     EFFZ=interp1(z,effectivedepth,Z,'nearest','extrap');
     EFFZprev=interp1(z,effectivedepth,Zprev,'nearest','extrap');
     DEN=(EFFZprev-EFFZ)./(Zprev-Z);
     tstep=time(n+1)-time(n);
     Cstep=zeros(length(zs),length(erosion));
     for ni=1:length(P)
         Cstep=Cstep+P(ni)./(l+(erosion.*DEN)./L(ni)).*...
             exp(-EFFZ./L(ni)).*...
             (1-exp(-tstep.*(l+(erosion.*DEN)./L(ni))));
     end
     C=C+Cstep.*exp(-l*time(n));
 end
 
 
end

