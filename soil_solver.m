function [  ] = soil_solver( profilefile )


disp(['-------------------------------------------------------------------------------------------------'])
disp(['CoSOILcal:' char(9) 'Cosmogenic Erosion Rate Calculator for samples in a variable-density profile v1.0'])
disp(['Written by:' char(9) 'Angel Rodes (angelrodes-at-gmail.com), CIAF, SUERC, Scotland. June 2019'])
disp(['File: ' profilefile])
disp(['Time: ' datestr(now)])



%% Input data

ndatasamples=csvread(profilefile,1,0);
ndatasamples(ndatasamples==0)=NaN;

numlines=size(ndatasamples,1);
ndens=0;
nconc=0;
for sample=1:1:numlines
    if ~isnan(ndatasamples(sample,2))
        ndens=ndens+1;
        depthd(ndens)=ndatasamples(sample,1);
        density(ndens)=ndatasamples(sample,2);
    end
    if ~isnan(ndatasamples(sample,3))
        nconc=nconc+1;
        sampledepth(nconc)=ndatasamples(sample,1);
        C(nconc)=ndatasamples(sample,3);
        dC(nconc)=ndatasamples(sample,4);
        nuclide(nconc)=ndatasamples(sample,5);
    end
end

site_lat=ndatasamples(1,6);
site_lon=ndatasamples(1,7);
site_elv=ndatasamples(1,8);
shielding=ndatasamples(1,9);
age=ndatasamples(1,10);

% display concentration data
disp(' ')
disp('Concentrations:')
for n=1:length(sampledepth)
disp([num2str(sampledepth(n)) ' cm' char(9) num2str(C(n)) ' +/- ' num2str(dC(n)) ' atoms(' num2str(nuclide(n)) ')/g'])
end

%% Make profile
z=unique(sort([0,depthd,sampledepth,logspace(0,log10(max(100000,10*max([depthd,sampledepth]))),100)]));
densityz=interp1(depthd,density,z,'nearest','extrap');
effectivedepth=cumsum(diff([0,z]).*densityz);

% plot density profile
figure
subplot(1,3,1)
hold on
plot(densityz,-z,'-k')
plot(density,-depthd,'.r')
ylim([-1.2*max([depthd,sampledepth]) 0])
xlim([0 1.2*max(densityz)])
set(gca, 'XAxisLocation', 'top')
xlabel('Density (g/cm3)')
ylabel('Depth (cm)')
% title('a')

% plot concentrations
subplot(1,3,2)
hold on
for n=1:length(sampledepth)
    plot(C(n),-sampledepth(n),'.r')
    plot([C(n)-dC(n),C(n)+dC(n)],[-sampledepth(n),-sampledepth(n)],'-r')
end
ylim([-1.2*max([depthd,sampledepth]) 0])
xlim([0 1.2*max(C+dC)])
set(gca, 'XAxisLocation', 'top')
xlabel('Concentrations (atoms/g)')
% ylabel('Depth (cm)')

%% Get production rates
for n=1:length(sampledepth)
Pr=Production_rate(site_lat,site_lon,site_elv,shielding,nuclide(n));
P(:,n)=Pr.P';
L(:,n)=Pr.L';
l(n)=Pr.l;
end

% display surface production rates
disp(' ')
disp('Surface production rates:')
for nuc=unique(nuclide)
Pr=Production_rate(site_lat,site_lon,site_elv,shielding,nuc);
disp([Pr.label char(9) num2str(sum(Pr.P)) ' atoms/g/year'])
end

%% Fit model to data
nresolution=30;
niter=10;

er=logspace(-6,1,nresolution);

convergence=0;
iter=0;
while convergence==0
    iter=iter+1;
    dev=desvmodel(P,L,l,z,effectivedepth,age,er,sampledepth,C,dC);
    ersol=[];
    for sample=1:size(dev,1)
        [~,selunique,~]=unique(dev(sample,:));
        ersol=[ersol,interp1(dev(sample,selunique),er(selunique),[-10,0,10],'linear','extrap')];
    end
    if length(unique(ersol))>2 && iter<niter
        er=logspace(log10(max(1e-6,min(ersol))),log10(max(ersol)),nresolution);
    else
        convergence=1;
    end
end
nresolution=1000;
er=logspace(log10(max(1e-6,min(ersol))),log10(max(ersol)),nresolution);
dev=desvmodel(P,L,l,z,effectivedepth,age,er,sampledepth,C,dC);

chisq=sum(dev.^2,1);

best=find(chisq==min(chisq),1,'first');
selonesigma=chisq<min(chisq)+length(C);
if sum(selonesigma)<3
    selonesigma=(min(abs(dev),1)<=1);
end

% display results
disp(' ')
disp(['Results for given age of ' num2str(age/1000) ' ka:'])
disp(['Best fit:' char(9) num2str(er(best)*1e4) ' m/Ma'])
disp(['One sigma resutls:' char(9) num2str(min(er(selonesigma))*1e4) ' - ' num2str(max(er(selonesigma))*1e4) ' m/Ma'])
disp(['Chi-squared values:' char(9) num2str(min(chisq(selonesigma))) ' - ' num2str(max(chisq(selonesigma))) ' '])

averresult=er(best);
uncertapprox=(max(er(selonesigma))-min(er(selonesigma)))/2;
precis=-(floor(log10(uncertapprox))-1);
averresult=round(averresult*10^precis)*10^(-precis);
uncertapprox=round(uncertapprox*10^precis)*10^(-precis);

%% plot results
subplot(1,3,2)
hold on

for n=1:length(C)
    Cmodel=model(P(:,n),L(:,n),l(n),z,effectivedepth,age,er(best),z);
    plot(Cmodel,-z,'-k')
    Cmodel=model(P(:,n),L(:,n),l(n),z,effectivedepth,age,max(er(selonesigma)),z);
    plot(Cmodel,-z,'--k')
    Cmodel=model(P(:,n),L(:,n),l(n),z,effectivedepth,age,min(er(selonesigma)),z);
    plot(Cmodel,-z,'--k')
end

for n=1:length(sampledepth)
    plot(C(n),-sampledepth(n),'.r')
    plot([C(n)-dC(n),C(n)+dC(n)],[-sampledepth(n),-sampledepth(n)],'-r')
end
text(min(C),-1.1*max(sampledepth),[num2str(averresult*1e4) '+/-' num2str(uncertapprox*1e4) ' m/Ma'])
ylim([-1.2*max([depthd,sampledepth]) 0])
xlim([0 1.2*max(max(C+dC),max(Cmodel)/1.2)])
set(gca, 'XAxisLocation', 'top')
xlabel('Concentrations (atoms/g)')
% title('b')

%% plot probability
subplot(1,3,3)
hold on
prob=exp(-chisq/2);
plot(er*1e4,prob,'-k')
text(averresult*1e4,max(prob)*1.05,['min chi-sq.=' num2str(min(round(chisq*100)/100))])
ylim([0 max(prob)*1.15])
xlim([(averresult-uncertapprox*5) (averresult+uncertapprox*5)]*1e4)
set(gca, 'YTick', []);
box on
xlabel('Erosion rate PDD (m/Ma)')
% title('c')


end

