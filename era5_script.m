%system('dir *.nc > yr.txt');
clear all;clc;
fid=fopen('yr1.txt');
C=textscan(fid, '%s')
fnames1 = C{1,1};

lat=double(ncread('era5.pl.monthly.1956.nc','latitude'));
lon=double(ncread('era5.pl.monthly.1956.nc','longitude'));
[lat,lon]=meshgrid(lat,lon);
lon1=lon(:,31:241);lat1=lat(:,31:241);    
%lat 101 111, lon 121 141
for m=1:72
new_name = [fnames1{m,1}];
u11=double(ncread(new_name,'t')); %4 850 7 700, 11 500, 17 200;
%u12=squeeze(u11(:,:,7,:));
%if size(u12,3)==730
%u13=squeeze(u12(:,31:241,546+1:546+122));
%else u13=squeeze(u12(:,31:241,548+1:548+122));
%end
%u14=squeeze(nanmean(reshape(u13,[720,211,2,61]),3));
%500 300 200; %4 850 7 700
%u15(:,:,:,m)=u14;
u12=squeeze(u11(:,31:241,1:2:17,:));
%u13=squeeze(nanmean(reshape(u12,[720,211,9,12]),3));
%u13=squeeze(reshape(u12,[720,361,3,12]));
%u14(:,:,:,:,m)=u13;
u13(:,:,:,:,m)=u12;
count =m
end

for l=1:60
for i=1:720
    for j=1:211
        tmp21=squeeze(u14(i,j,l,:));
        tmp22=detrend(tmp21,1);
        olrjf_dt(i,j,l,:)=tmp22;
    end
end
count=l
end


for m=1:72
new_name = [fnames1{m,1}];
u11=double(ncread(new_name,'z'))/9.81; %4 850 7 700, 11 500, 17 200;
u12=squeeze(u11(121:141,101:111,17,:));
u13(:,:,:,m)=u12;
count =m
end

%pr2=reshape(pr1,[720,361,12,84]);
%pr21=pr2(:,:,:,10+1:10+72);
z11=u13(:,:,:,:,31)

for l=1:12
for k=1:9
for i=1:720
    for j=1:211
        tmp21=squeeze(u13(i,j,k,l,:));
        tmp22=detrend(tmp21,1);
        t_dt(i,j,k,l,:)=tmp22;
    end
end
end
count=l
end

for l=1:12
for i=1:720
    for j=1:211
        tmp21=squeeze(olr12(i,j,l,:));
        tmp22=detrend(tmp21,1);
        olr_dt(i,j,l,:)=tmp22;
    end
end
count=l
end


ts1=nino34_cpc';
for i=1:12
    tmp31=squeeze(ts1(i,:));
    tmp32=detrend(tmp31,1);
    ts_dt(i,:)=tmp32;
end

for l=1:12
for k=1:5
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(9,31+1:31+40), squeeze(z_dt(i,j,k,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k,l)=r(1,2);
        else corr_sst(i,j,k,l)=nan;
        end
    end
end
end
count=l
end

month1={'January','February','March','April','May', 'June', 'July', 'August', 'September', 'October', 'November' ,'December'};

%close all;
for i=1:5
figure();h=pcolor(lon1,lat1,corr_sst(:,:,3,i+7));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,3,i+7),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-Z500, lag-0, ', num2str(month1{i+7})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=1.5*500*(30/30);height=1.5*300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*5*(30/30) 1.5*3 ]);
%print('-dpng', '-r400', ['corr_ninoz500_lag0_',num2str(i)]);
count=i
end

close all;
figure();h=pcolor(lon1,lat1,u14_jfm);shading flat;colorbar;colormap(jet(30));caxis([-60 60]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,u14_jfm,[40 60 80],'color','k','LineWidth',1);
title(['U200 (ms^{-1}), JFM'],'FontSize',20);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=1.5*300*(150/90);height=1.5*300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*3*(150/90) 1.5*3 ]);
print('-dpng', '-r400', ['u_jfm']);
count=i

%tgrad  0-50 (0 101) -5  5 (111 131)
% tgrad 0-50 (0 101) 35 45 ( 31 - 51)

for l=1:12
for k=1:3
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(10,31+1-1:31+40-1), squeeze(z_dt(i,j,k,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k,l)=r(1,2);
        else corr_sst(i,j,k,l)=nan;
        end
    end
end
end
count=l
end

close all;
for i=1:3
figure();h=pcolor(lon1,lat1,corr_sst(:,:,3,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%ylim([-60 60]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,3,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-V200 lag: Month ', num2str(i)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=600*(40/30);height=300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 5 ]);
print('-dpng', '-r400', ['corr_prnao_lag5_ndjfm_',num2str(i)]);
count=i
end

z_dt_2=squeeze(mean(z_dt,3));

for l=1:12
%for k=1:3
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(10,31+1:31+40), squeeze(z_dt_2(i,j,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
%end
count=l
end

%for l=1:12
%for k=1:3
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(10,31+1:31+40), squeeze(rv_nv(i,j,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end
%end
%count=l
%end

close all;
for i=1:12
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-T200, lag-3, ', num2str(month1{i})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=1.5*500*(30/30);height=1.5*300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*5*(30/30) 1.5*3 ]);
%print('-dpng', '-r400', ['corr_ninot_lag3_',num2str(i)]);
count=i
end
figure();h=pcolor(lon1,lat1,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst,[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-T200, lag-3, ', num2str(month1{10})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=1.5*500*(30/30);height=1.5*300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*5*(30/30) 1.5*3 ]);
%print('-dpng', '-r400', ['corr_ninot_lag3_',num2str(i)]);

%vort 30-35N, 40-60E
%u200 40-60E, 25-30N
%u200 40-60E, 30-35N

close all;
for i=1:3
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%ylim([-60 60]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-V200 lag: Month ', num2str(i)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=600*(40/30);height=300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 5 ]);
%print('-dpng', '-r400', ['corr_prnao_lag5_ndjfm_',num2str(i)]);
count=i
end



for m=1:72
new_name = [fnames1{m,1}];
u11=double(ncread(new_name,'u')); %4 850 7 700, 11 500, 17 200;
%u12=squeeze(u11(:,61:241,[7,11,17],:));
%u13(:,:,:,:,m)=u12;
u12=squeeze(u11(:,61:241,[4,17],:));
u121=squeeze(u12(:,:,2,:)-u12(:,:,1,:));
u13(:,:,:,m)=u121;
count =m
end

for l=1:12
for i=1:720
    for j=1:181
        tmp21=squeeze(u13(i,j,l,:));
        tmp22=detrend(tmp21,1);
        z_dt(i,j,l,:)=tmp22;
    end
end
count=l
end

ts1=nino34_cpc';
for i=1:12
    tmp31=squeeze(ts1(i,:));
    tmp32=detrend(tmp31,1);
    ts_dt(i,:)=tmp32;
end

for l=1:12
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(l,31+1:31+40), squeeze(z_dt(i,j,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
count=l
end


month1={'January','February','March','April','May', 'June', 'July', 'August', 'September', 'October', 'November' ,'December'};

for l=1:3
%for k=1:3
for i=1:720
    for j=1:181
        [r,p] = corrcoef(ts_dt(12,31+1-1:31+40-1), squeeze(z_dt(i,j,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
%end
count=l
end

close all;
for i=1:12
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-T200, lag-0, ', num2str(month1{i})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=500*(30/30);height=300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(30/30) 3 ]);
%print('-dpng', '-r400', ['corr_ninot_lag0_',num2str(i)]);
count=i
end


rv1=z_dt_2(81:121,51:61,:,:);

rv2=squeeze(mean(mean(rv1,1),2));

for l=1:12
%for k=1:3
for i=1:41
    for j=1:41
        [r,p] = corrcoef(rv2(l,31+1:31+40), squeeze(pr_hma_dt(i,j,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
%end
count=l
end

close all;
for i=1:12
figure();h=pcolor(lon_hma,lat_hma,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst(:,:,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['HMA PR-Vort, lag-0, ', num2str(month1{i})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=500*(40/30);height=500; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 4 ]);
print('-dpng', '-r400', ['corr_rvhma_lag0_',num2str(i)]);
count=i
end


rv1=z_dt(111:131,21:31,:,:,:); 

%6070E 121-141
%30-35 51-61

%30-35 35-40,60-70
%45-50, 55-65 60-70

rv2=squeeze(mean(mean(rv1,1),2));

for l=1:12
%for k=1:3
for i=1:41
    for j=1:41
        [r,p] = corrcoef(rv2(2,l,31+1:31+40), squeeze(pr_hma_dt(i,j,l,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
%end
count=l
end

close all;
for i=1:2
figure();h=pcolor(lon_hma,lat_hma,corr_sst(:,:,i+10));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst(:,:,i+10),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['HMA PR-Z500, lag-0, ', num2str(month1{i+10})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=500*(40/30);height=500; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 4 ]);
print('-dpng', '-r400', ['corr_z500hma2_lag0_',num2str(i)]);
count=i
end

%%corr era5precip
pr2=reshape(pr1,[720,361,12,84]);
pr21=pr2(:,:,:,10+1:10+72);
pr22=pr21(:,61:241,:,:);

for l=1:12
for i=1:720
    for j=1:181
        tmp21=squeeze(pr22(i,j,l,:));
        tmp22=detrend(tmp21,1);
        z_dt(i,j,l,:)=tmp22;
    end
 end
count=l
end

ts1=nino34_cpc';
for i=1:12
    tmp31=squeeze(ts1(i,:));
    tmp32=detrend(tmp31,1);
    ts_dt(i,:)=tmp32;
end

lon2=lon1(123:163,33:73);lat2=lat1(123:163,33:73);
z_dt2=z_dt(123:163,33:73,:,:);

pr_nd_era_dt=squeeze(z_dt2(:,:,11,:)+z_dt2(:,:,12,:))*1000;

for l=1:12
for i=1:41
    for j=1:41
        [r,p] = corrcoef(ts2(:,1), squeeze(pr_nd_era_dt(i,j,31+1-10:31+40-10)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
count=l
end

month1={'January','February','March','April','May', 'June', 'July', 'August', 'September', 'October', 'November' ,'December'};

close all;
for i=1:6
figure();h=pcolor(lon2,lat2,corr_sst(:,:,i+6));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon2,lat2,corr_sst(:,:,i+6),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-HMAPR, ', num2str(month1{i+6})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=700*(40/30);height=700; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7*(40/30) 7 ]);
print('-dpng', '-r400', ['corr_hmaninoera_',num2str(i+6)]);
count=i
end

for l=1:12
for i=1:41
    for j=1:41
        [r,p] = corrcoef(ts2(:,1), squeeze(pr_nd_era_dt(i,j,31+1-10:31+40-10)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
count=l
end

month1={'January','February','March','April','May', 'June', 'July', 'August', 'September', 'October', 'November' ,'December'};

close all;
figure();h=pcolor(lon2,lat2,corr_sst(:,:,1));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon2,lat2,corr_sst(:,:,1),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-HMAPR, Nov+Dec'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=700*(40/30);height=700; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7*(40/30) 7 ]);
print('-dpng', '-r400', ['corr_hmaninond_era_nd']);
count=i
%end


%
%partial correlation

%pc(pr(t),enso(t-1)|tgrad(t-1),stj(t),vort1(t),vort2(t),twio(t-1),ivt1(t),ivt2(t),z200(t),slp(t))
ts1=nino3';
for i=1:12
    tmp31=squeeze(ts1(i,:));
    tmp32=detrend(tmp31,1);
    ts_dt(i,:)=tmp32;
end

ts1=nino34_cpc';
ts2=ts1(:,31+1:31+40);
for i=1:12
    tmp31=squeeze(ts2(i,:));
    tmp32=detrend(tmp31,1);
    ts_dt(i,:)=tmp32;
end

pr_hma1=pr_hma(:,:,:,31+1:31+40);
for k=1:12
    for i=1:41
        for j=1:41
        tmp21=squeeze(pr_hma1(i,j,k,:));
        tmp22=detrend(tmp21,1);
        pr_hma_dt(i,j,k,:)=tmp22;
        end
    end
    count=k
end
for k=1:12
    for i=1:41
        for j=1:41
        tmp21=squeeze(pr_hma(i,j,k,:));
        tmp22=detrend(tmp21,1);
        pr_hma_dt(i,j,k,:)=tmp22;
        end
    end
    count=k
end

for l=1:12
%for k=1:3
for i=1:41
    for j=1:41
        [r,p] = corrcoef(ts_dt(l,31+1:31+40), squeeze(pr_hma_dt(i,j,l+3,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
%end
count=l
end

for l=1:12
for i=1:41
    for j=1:41
        [r,p] = corrcoef(ts_dt(l,31+1:31+40), squeeze(pr_ond_dt(i,j,31+1:31+40)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
end

for l=1:12
for i=1:41
    for j=1:41
        [r,p] = corrcoef(ts_dt(l,:), squeeze(pr_hma_dt(i,j,l+1,:)));
        if p(1,2)<=0.05;
            corr_sst(i,j,l)=r(1,2);
        else corr_sst(i,j,l)=nan;
        end
    end
end
count=l
end

month1={'January','February','March','April','May', 'June', 'July', 'August', 'September', 'October', 'November' ,'December'};

%jfm precip nino3
%ond precip

close all;
for i=1:12
figure();h=pcolor(lon_hma,lat_hma,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);%xlim([0 150]);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst(:,:,i),[-0.9, -0.7, -0.5, 0.5, 0.7, 0.9],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ENSO-Z500, lag-0, ', num2str(month1{i})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=100;width=1.5*500*(30/30);height=1.5*300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 1.5*5*(30/30) 1.5*3 ]);
%print('-dpng', '-r400', ['corr_ninoz500_lag0_',num2str(i)]);
count=i
end


rv_nv=mean(rv_2(11:12,:),1);
uu_nv=mean(u_ut(11:12,:),1);
ul_nv=mean(u_lt(11:12,:),1);


%lag=5;
for k=1:4
for i=1:41
    for j=1:41
        A=[ts_dt(10,1+10*(k-1)+1:40+10*(k-1)+1),squeeze(pr_nd_dt(i,j,1+10*(k-1)+1:40+10*(k-1)+1))];
        B=[dt1(10,1+10*(k-1)+1:40+10*(k-1)+1),uu_nv(1,1+10*(k-1)+1:40+10*(k-1)+1),ul_nv(1,1+10*(k-1)+1:40+10*(k-1)+1),rv_nv(1,1+10*(k-1)+1:40+10*(k-1)+1)];
        [r,p] = partialcorr(A,B);
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end

for i=1:41
    for j=1:41
        A=[ts_dt(31+1:31:40,10),squeeze(pr_nd_dt(i,j,31+1:31:40))];
%        B=[dt1(31+1:31:40,10),uu_nv(31+1:31:40,1),ul_nv(31+1:31:40,1),rv_nv(31+1:31:40,1)];
        B=[dt1(31+1:31:40,10)];
        [r,p] = partialcorr(A,B);
%        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
%        else corr_sst(i,j)=nan;
%        end
    end
end
close all;
figure();h=pcolor(lon_hma,lat_hma,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst,[-0.4, 0.4],'color','k','LineWidth',1);
%title(['IVT-HMAPR|IOD: Lag-5 ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=500*(40/30);height=500; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 5 ]);
%print('-dpng', '-r400', ['parcorr_2330_privt0_iod5_lag5_ndjfm_',num2str(i)]);


close all;
for i=1:4
figure();h=pcolor(lon_hma,lat_hma,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst(:,:,i),[-0.4, 0.4],'color','k','LineWidth',1);
title(['IVT-HMAPR|IOD: Lag-5 ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=500*(40/30);height=500; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5*(40/30) 5 ]);
%print('-dpng', '-r400', ['parcorr_2330_privt0_iod5_lag5_ndjfm_',num2str(i)]);
count=i
end

close all;
for i=1:42
figure();h=pcolor(lon_hma,lat_hma,pr_nd_dt(:,:,i+30));shading flat;colorbar;colormap(jet(40));caxis([-10 10]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
count=i
end

ivt121=squeeze(mean(ivt12(:,11:12),2));
ivt221=squeeze(mean(ivt22(:,11:12),2));
ivt321=squeeze(mean(ivt32(:,11:12),2));
twio1=squeeze(mean(twio(:,11:12),2));

%y1=-10;y2=-10;%1991-2020
%y1=-10+30;y2=-10;%1961-1990
y1=10;y2=0;%1970-2010;
%y1=20;y2=0;%1960-2000;

z11=mean(z_cgt(11:12,:),1)';

z1=   z11(31+1-y1:31+40-y1+y2,1);
t1=   dt1(31+1-y1:31+40-y1+y2,:);
ts1=ts_dt(31+1-y1:31+40-y1+y2,:);
u1= uu_nv(31+1-y1:31+40-y1+y2,1);
u2= ul_nv(31+1-y1:31+40-y1+y2,1);
rv1=rv_nv(31+1-y1:31+40-y1+y2,1);
ivt1=ivt121(31+1-y1:31+40-y1+y2,1);
ivt2=ivt221(31+1-y1:31+40-y1+y2,1); %(45-65, 23-30N)
ivt3=ivt321(31+1-y1:31+40-y1+y2,1);
twio1= twio(31+1-y1:31+40-y1+y2,:);
%u1= u_ut(1,31+1-1:31+40-1)';
%u2= u_lt(1,31+1-1:31+40-1)';
%rv1=rv_2(1,31+1-1:31+40-1);
ts2=mean(ts1(:,11:12),2);
t2=mean(t1(:,11:12),2);
twio2=mean(twio1(:,11:12),2);;

X = [t1(:,11),u1,u2,rv1,ivt1];%t1(:,11),u1,u2,rv1,ivt2,twio1(:,11), z1;
mdl_nino = fitlm(X,ts1(:,10));
e_nino=mdl_nino.Residuals.Raw;


warning("off");
for i=1:41
    for j=1:41
        pr1=squeeze(pr_nd_dt(i,j,31+1-y1:31+40-y1+y2));
        mdl_pr=fitlm(X,pr1);
        e_pr1=mdl_pr.Residuals.Raw;
        e_pr2(i,j,:)=e_pr1;
    end
end

for i=1:41
    for j=1:41
        [r,p] = corrcoef(e_nino(:,1),e_pr2(i,j,:));
        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end


%close all;
figure();h=pcolor(lon_hma,lat_hma,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst,[-0.4, 0.4],'color','k','LineWidth',1);
title(['ParCorr: [HMAPR-ENSO|IVT], 1970-2010'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=400*(40/30);height=400; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7*(40/30) 7 ]);
%print('-dpng', '-r400', ['parcorr_era5_prninond_ivt_lag0']);
