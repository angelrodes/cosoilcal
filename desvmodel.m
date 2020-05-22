function [ devs ] = desvmodel( P,L,l,z,effectivedepth,age,erosion,sampledepth,C,dC)
% Calculates deviations as (C-Model)/dC

% resize paramanters
devs=zeros(length(C),length(erosion))+NaN;
C=C(:);
dC=dC(:);


for n=1:length(C)
    devs(n,:)=(C(n)-model(P(:,n),L(:,n),l(n),z,effectivedepth,age,erosion,sampledepth(n)))./dC(n);
end

end

