clear all;clc;
fid=fopen('yr1.txt');
C=textscan(fid, '%s')
fnames1 = C{1,1};

lat=double(ncread('era5.pl.monthly.1940.nc','latitude'));
lon=double(ncread('era5.pl.monthly.1940.nc','longitude'));
[lat,lon]=meshgrid(lat,lon);
lon1=lon(:,31:241);lat1=lat(:,31:241);    
for m=1:72
new_name = [fnames1{m,1}];
u11=double(ncread(new_name,'u'));
u12=squeeze(u11(:,31:241,[4,17],:));%500 300 200;
u13(:,:,:,:,m)=u12;
count =m
end

for l=1:12
for k=1:2
for i=1:720
    for j=1:211
        tmp21=squeeze(u13(i,j,k,l,:));
        tmp22=detrend(tmp21,1);
        u_dt(i,j,k,l,:)=tmp22;
    end
end
end
count=l
end

sst1=sst3(:,:,:,50:50+71);
for l=1:12
for i=1:360
    for j=1:180
        tmp21=squeeze(sst1(i,j,l,:));
        tmp22=detrend(tmp21,1);
        sst_dt(i,j,l,:)=tmp22;
    end
end
count=l
end


for l=1:70
for i=1:41
for j=1:41
pr_hma_jfm(i,j,l)=pr_hma(i,j,1,l+1)+pr_hma(i,j,2,l+1)+pr_hma(i,j,3,l+1); %(1951-2020)
end
end
count=l
end

for l=1:70
for i=1:41
for j=1:41
pr_hma_nd(i,j,l)=pr_hma(i,j,11,l+1)+pr_hma(i,j,12,l+1); %(1951-2020)
end
end
count=l
end

for i=1:41
    for j=1:41
        tmp21=squeeze(pr_hma_jfm(i,j,:));
        tmp22=detrend(tmp21,1);
        pr_jfm_dt(i,j,:)=tmp22;
    end
end

for i=1:41
    for j=1:41
        tmp21=squeeze(pr_hma_nd(i,j,:));
        tmp22=detrend(tmp21,1);
        pr_nd_dt(i,j,:)=tmp22;
    end
end

%29.75-43.25N, 62.25-78.25E
%12-39,        3-35

%36.25N 25

%34-43N 21-39, 65-75E 9 29;


pr11=pr_jfm_dt(3:35,12:25,:);
pr12=squeeze(nanmean(nanmean(pr11,1),2));
pr1_jfm=pr12;

pr21=pr_nd_dt(3:35,12:39,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_nd=pr22;

pr21=pr_nd_dt(3:35,25:39,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_nd_n=pr22;

pr21=pr_nd_dt(3:35,12:25,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_nd_s=pr22;

%cat
for i=1:720
    for j=1:211
        if lon1(i,j)>180;
        lon11(i,j)=lon1(i,j)-360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:361,:);lon13=lon11(362:720,:);
lon14=cat(1,lon13,lon12);lat14=lat1;

yr=40;
for k=1:4
for i=1:720
    for j=1:211
        [r,p] = corrcoef(pr1_nd(1+10*(k-1):yr+10*(k-1),1), squeeze(z_dt(i,j,3,11,1+10*(k-1)+1:yr+10*(k-1)+1)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end
%corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for i=1:4
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,corcoastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.4, 0.4],'color','k','LineWidth',1);
title(['HMAPR-Z500 lag-0: ND ', num2str(1951+(i-1)*10),'-',num2str(1950+yr+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=200*(360/105);height=200; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_prz500_',num2str(yr),'_lag0_nd_',num2str(i)]);
count=i
end


for i=1:360
    for j=1:180
        if lon1(i,j)<0;
        lon11(i,j)=lon1(i,j)+360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lon14=cat(1,lon13,lon12);lat14=lat1;

sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11); 


yr=40;
for k=1:4
for i=1:360
    for j=1:180
        [r,p] = corrcoef(pr1_nd(1+10*(k-1):yr+10*(k-1),1), squeeze(sst_dt(i,j,10,1+10*(k-1)+1:yr+10*(k-1)+1)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end
corr_sst11=corr_sst(1:180,:,:);corr_sst12=corr_sst(181:360,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for i=1:4
figure();h=pcolor(lon14,lat14,corr_sst2(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,corr_sst2(:,:,i),[-0.4, 0.4],'color','k','LineWidth',1);
title(['HMAPR-SST lag-0: ND ', num2str(1951+(i-1)*10),'-',num2str(1950+yr+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;ylim([-30 60]);
x0=100;y0=50;width=300*(360/95);height=300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/95) 3 ]);
%print('-dpng', '-r400', ['corr_prsst_',num2str(yr),'_lag0_nd_',num2str(i)]);
count=i
end


lon1=double(ncread('data_1.nc','longitude'));
lat1=double(ncread('data_1.nc','latitude'));
[lat1,lon1]=meshgrid(lat1,lon1);
p1=double(ncread('data_1.nc','msl'));
p2=reshape(p1,[720,361,12,84]);

p21=p2(:,:,:,11:11+71);
p22=p21(:,31:241,:,:);

for k=1:12
for i=1:720
for j=1:211
tmp21=squeeze(p22(i,j,k,:));
tmp22=detrend(tmp21,1);
slp_dt(i,j,k,:)=tmp22;
end
end
count=k
end

lon1=double(ncread('data_2.nc','longitude'));
lat1=double(ncread('data_2.nc','latitude'));
[lat1,lon1]=meshgrid(lat1,lon1);
p1=double(ncread('data_2.nc','ttr'));
p2=reshape(p1,[720,361,12,84]);

p21=p2(:,:,:,11:11+71);
p22=p21(:,31:241,:,:);
for k=1:12
for i=1:720
for j=1:211
tmp21=squeeze(p22(i,j,k,:));
tmp22=detrend(tmp21,1);
olr_dt(i,j,k,:)=tmp22;
end
end
count=k
end

for i=1:720
    for j=1:211
        if lon1(i,j)>180;
        lon11(i,j)=lon1(i,j)-360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:361,:);lon13=lon11(362:720,:);
lon14=cat(1,lon13,lon12);lat14=lat1;



for i=1:720
    for j=1:361
        [r,p] = corrcoef(pr2(11,2:71), squeeze(zl875_dt(i,j,1,10,2:71)));
        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end
%corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 

%close all;
figure();h=pcolor(lon,lat,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon,lat,corr_sst,[-0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-vimd'],'FontSize',15);ylim([-45 45]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_n_prz200_cov_pc2']);

close all;
for k=1:12
figure();h=pcolor(lon14,lat14,corr_sst2(:,:,k));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,corr_sst2(:,:,k),[-0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-Z200: ',month1{k}],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_n_prz200_cov_pc2']);
end

for i=1:360
    for j=1:180
        [r,p] = corrcoef(pc(1,:), squeeze(sst_dt(i,j,10,2:71)));
        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end
%corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for k=1:12
figure();h=pcolor(lon1,lat1,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst,[-0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-SST'],'FontSize',15);ylim([-30 75]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_n_prsst_cov_pc2']);
end


for i=1:720
    for j=1:211
        [r,p] = corrcoef(pc(3,:), squeeze(u_dt(i,j,4,3,2:71)));
        if p(1,2)<=0.05;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end
%corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
figure();h=pcolor(lon1,lat1,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst,[-0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-Z200: '],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_n_prz200_cov_pc2']);


close all;
for i=1:8
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.4, 0.4],'color','k','LineWidth',1);
title(['HMAPR-Z500 lag-0: ND ', num2str(1951+(i-1)*5),'-',num2str(1950+yr+(i-1)*5)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=200*(360/105);height=200; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_prz500_',num2str(yr),'_lag0_nd_',num2str(i)]);
count=i
end



%1971-2010
%z1 enso: 180 240, -20-20:  361-481, 111 191
%z2 275-335, 26-46:         551-671, 59-99 % 561-661, 59 99
%z3 200-340, 50-75:         401-681, 1-51; 1 % 511 681
%z3 250-340, 50-75:         501-681, 1-51; 2
%z3 300-340, 50-75:         601-681, 1-51; 3
%z4 2:30, 37-67:            5:61, 17:77
%z5 35-70, 25-50:           71:141,  51:101          
%z6 80:115, 18:36:          161:231, 79:115
%z7: 120-145, 35:50:        241:291, 51:81 
%sst1: 165:235, -15:15:     166:236, 75:106
z11=squeeze(z_dt(:,:,3,11,:));

corr31=corr_sst(:,:,3);
for i=1:720
    for j=1:211
        corr32(i,j)=corr31(i,j)*0;
    end
end

for i=1:72
z12(:,:,i)=z11(:,:,i)+corr32(:,:);
end
%z1 enso: 180 240, -20-20:  361-481, 111 191
%z2 275-335, 26-46:         551-671, 59-99 % 561-661, 59 99
%z3 200-340, 50-75:         401-681, 1-51; 1 % 511 681
%z3 250-340, 50-75:         501-681, 1-51; 2
%z3 300-340, 50-75:         601-681, 1-51; 3
%z4 2:30, 37-67:            5:61, 17:77
%z5 35-70, 25-50:           71:141,  51:101          
%z6 80:115, 18:36:          161:231, 79:115
%z7: 120-145, 35:50:        241:291, 51:81 
%sst1: 170:230, -10:10:     171:231, 80:101

z1 =squeeze(nanmean(nanmean(z12(361:481, 111:191,:),1),2));
z2 =squeeze(nanmean(nanmean(z12(551:671, 59:99,:),1),2));
z3 =squeeze(nanmean(nanmean(z12(601:681, 1:51,:),1),2));
z4 =squeeze(nanmean(nanmean(z12(5:61, 17:77,:),1),2));
z5 =squeeze(nanmean(nanmean(z12(71:141,  51:101,:),1),2));          
z6 =squeeze(nanmean(nanmean(z12(161:231, 79:115,:),1),2));
z7 =squeeze(nanmean(nanmean(z12(241:291, 51:81,:),1),2));
 

z31=z11(401:681, 1:51,:);
latz31=lat1(401:681, 1:51,:);
z32=z31.cosd(latz31);
%z33=(squeeze(sum(sum(z32,1),2)))/(squeeze(sum(sum(cosd(latz31),1),2)));
z33=(squeeze(mean(mean(z32,1),2)));

for i=1:360
    for j=1:180
        if lon1(i,j)<0;
        lon11(i,j)=lon1(i,j)+360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lon14=cat(1,lon13,lon12);lat14=lat1;
sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11);

sst21=squeeze(sst2(171:231, 80:101,10,:));

lont1=lon14(166:236,75:106);
latt1=lat14(166:236,75:106);

sst1=squeeze((mean(mean(sst21,1),2)));% (october)

%corrcoef
[r1,p1]=corrcoef(pr1_nd(21:60,1),sst1(22:61,1));
rsst=r1(1,2);psst=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z1(22:61,1));
rz1=r1(1,2);pz1=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z2(22:61,1));
rz2=r1(1,2);pz2=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z3(22:61,1));
rz3=r1(1,2);pz3=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z4(22:61,1));
rz4=r1(1,2);pz4=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z5(22:61,1));
rz5=r1(1,2);pz5=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z6(22:61,1));
rz6=r1(1,2);pz6=p1(1,2);
[r1,p1]=corrcoef(pr1_nd(21:60,1),z7(22:61,1));
rz7=r1(1,2);pz7=p1(1,2);

%partial corr 
%pr-sst|z
[r1,p1]=partialcorr([pr1_nd,sst1(2:71,1)],[z1(2:71,1)])




%1951-2020
z11=squeeze(z_dt(:,:,3,11,:));

corr31=corr_sst;
for i=1:720
    for j=1:211
        corr32(i,j)=corr31(i,j)*0;
    end
end

for i=1:72
z12(:,:,i)=z11(:,:,i)+corr32(:,:);
end
%z1 enso: 170 250, -20-20:  341-501, 111 191
%z1 enso: 160 260, -25-25:  321-521, 101 201
%z2 270-330, 24-44:         541-661, 63-103
%z3 280-340, 50-75:         561-701, 1-51;
%z4 0:35, 40-75:            1:71, 1:71
%z5 35-75, 27-55:           71:151,  41:97          
%z6 86:113, 20:33:          173:227, 85:111
%z7: 115-145, 36:52:        231:291, 47:79
%sst1: 170:230, -10:15:     171:231, 75:101
%sst1: 165:235, -15:15:     166:236, 75:106
z1 =squeeze(nanmean(nanmean(z12(341:501, 111:191,:),1),2));
z2 =squeeze(nanmean(nanmean(z12(541:661, 63:103,:),1),2));
z3 =squeeze(nanmean(nanmean(z12(561:681, 1:51,:),1),2));
z4 =squeeze(nanmean(nanmean(z12(1:71, 1:71,:),1),2));
z5 =squeeze(nanmean(nanmean(z12(71:151,  41:97,:),1),2));          
z6 =squeeze(nanmean(nanmean(z12(173:227, 85:111,:),1),2));
z7 =squeeze(nanmean(nanmean(z12(231:291, 47:79,:),1),2));
 

%z31=z11(401:681, 1:51,:);
%latz31=lat1(401:681, 1:51,:);
%z32=z31.cosd(latz31);
%z33=(squeeze(sum(sum(z32,1),2)))/(squeeze(sum(sum(cosd(latz31),1),2)));
%z33=(squeeze(mean(mean(z32,1),2)));

for i=1:360
    for j=1:180
        if lon1(i,j)<0;
        lon11(i,j)=lon1(i,j)+360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lon14=cat(1,lon13,lon12);lat14=lat1;
sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11);
sst21=squeeze(sst2(:,:,10,:));

corr31=corr_sst2;
for i=1:360
    for j=1:180
        corr32(i,j)=corr31(i,j)*0;
    end
end

for i=1:72
sst22(:,:,i)=sst21(:,:,i)+corr32(:,:);
end


sst23=squeeze(sst22(171:231, 75:101,:));

lont1=lon14(171:241,70:101);
latt1=lat14(171:241,70:101);

sst1=squeeze((nanmean(nanmean(sst23,1),2)));% (october)

%corrcoef
[r1,p1]=corrcoef(pr1_nd,sst1(2:71,1));
rsst=r1(1,2);psst=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z1(2:71,1));
rz1=r1(1,2);pz1=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z2(2:71,1));
rz2=r1(1,2);pz2=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z3(2:71,1));
rz3=r1(1,2);pz3=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z4(2:71,1));
rz4=r1(1,2);pz4=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z5(2:71,1));
rz5=r1(1,2);pz5=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z6(2:71,1));
rz6=r1(1,2);pz6=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z7(2:71,1));
rz7=r1(1,2);pz7=p1(1,2);

%partial corr 
%pr-sst|z
[r1,p1]=partialcorr([pr1_nd,sst1(2:71,1)],[z1(2:71,1)])



%
close all;
for i=1:72
figure();h=pcolor(lon1,lat1,sst_dt(:,:,10,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon1,lat1,sst_dt(:,:,10,i),[-0.5, 0.5],'color','k','LineWidth',1);
title(['SSTA, Oct ', num2str(1949+i)],'FontSize',20);ylim([-30 45]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=50;y0=200;width=200*(300/80);height=200; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_prz500_',num2str(yr),'_lag0_nd_',num2str(i)]);
count=i
end


%
%enso years
k11=[1951,1953,1954,1955,1957,1963,1964,1965,1968,1969,1970,1971, 1972,1973,1974,1975,1976,1977,1979,1982,1983,1984,1986,1987,1988, 1991,1994,1995,1997,1998,1999,2000,2002,2004,2006,2007,2009,2010, 2011,2014,2015,2016,2017,2018,2020];
k12=k11-1950;
k21=[1951,1953,1957,1963,1965,1968,1969, 1972,1976,1977,1979,1982,1986,1987, 1991,1994,1997,2002,2004,2006,2009,2014,2015,2018];
k22=k21-1950;
k31=[1954,1955,1964,1970,1971,1973,1974,1975,1983,1984,1988,1995,1998,1999,2000,2007,2010, 2011,2016,2017,2020];
k32=k31-1950;

for i=1:24
    z_dt_en(:,:,:,:,i)=z_dt(:,:,:,:,k22(i)+1);
end
for i=1:21
    z_dt_ln(:,:,:,:,i)=z_dt(:,:,:,:,k32(i)+1);
end
for i=1:45
    z_dt_enso(:,:,:,:,i)=z_dt(:,:,:,:,k12(i)+1);
end

for i=1:24
    s_dt_en(:,:,:,i)=slp_dt(:,:,:,k22(i)+1);
end
for i=1:21
    s_dt_ln(:,:,:,i)=slp_dt(:,:,:,k32(i)+1);
end
for i=1:45
    s_dt_enso(:,:,:,i)=slp_dt(:,:,:,k12(i)+1);
end

for i=1:45
    pr1_nd_enso(i,1)=pr1_nd(k12(i),1);
end
for i=1:24
    pr1_nd_en(i,1)=pr1_nd(k22(i),1);
end
for i=1:21
    pr1_nd_ln(i,1)=pr1_nd(k32(i),1);
end

for i=1:720
    for j=1:211
        if lon1(i,j)>180;
        lon11(i,j)=lon1(i,j)-360;
        else lon11(i,j)=lon1(i,j);
        end
    end
end
lon12=lon11(1:361,:);lon13=lon11(362:720,:);
lon14=cat(1,lon13,lon12);lat14=lat1;

zdt11=z_dt_en(1:361,:,:,:,:);zdt12=z_dt_en(362:720,:,:,:,:);
z_dt_en_2=cat(1,zdt12,zdt11);
udt11=u_dt_en(1:361,:,:,:,:);udt12=u_dt_en(362:720,:,:,:,:);
u_dt_en_2=cat(1,udt12,udt11);
vdt11=v_dt_en(1:361,:,:,:,:);vdt12=v_dt_en(362:720,:,:,:,:);
v_dt_en_2=cat(1,vdt12,vdt11);
zdt11=z_dt_ln(1:361,:,:,:,:);zdt12=z_dt_ln(362:720,:,:,:,:);
z_dt_ln_2=cat(1,zdt12,zdt11);
udt11=u_dt_ln(1:361,:,:,:,:);udt12=u_dt_ln(362:720,:,:,:,:);
u_dt_ln_2=cat(1,udt12,udt11);
vdt11=v_dt_ln(1:361,:,:,:,:);vdt12=v_dt_ln(362:720,:,:,:,:);
v_dt_ln_2=cat(1,vdt12,vdt11);


long1=-180:4:180;latg1=-30:4:75;
[latg1,long1]=meshgrid(latg1,long1);
z_dt_2=squeeze(mean(z_dt_ln_2(:,:,2,11,:),5));%-squeeze(mean(z_dt_ln_2(:,:,1,11,5:17),5));
u_dt_2=squeeze(mean(u_dt_ln_2(:,:,2,11,:),5));%-squeeze(mean(u_dt_ln_2(:,:,1,11,5:17),5));
v_dt_2=squeeze(mean(v_dt_ln_2(:,:,2,11,:),5));%-squeeze(mean(v_dt_ln_2(:,:,1,11,5:17),5));

ugrid=griddata(lon14,lat14,u_dt_2,long1,latg1);
vgrid=griddata(lon14,lat14,v_dt_2,long1,latg1);

close all;figure();h=pcolor(lon14,lat14,z_dt_2);shading flat;
hold on; quiver(long1,latg1,ugrid,vgrid,'color','k','linewidth',1);
colormap(jet(40));caxis([-50 50]);ylim([-15 75]);
cbh=colorbar;set(cbh,'YTick',[-50:10:50]);
%cbh=colorbar;set(cbh,'YTick',[-25:10:25]);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,z_dt_2,[-20 20],'color','k','LineWidth',1);
title(['Z300, La Nina: Nov [1951-2020]'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
set (gca,'YTick',[-30 -15 0 15 30 45 60 75 90]);
set(gca,'XTick',[-60 0 60 120 180 240 300 360]);
x0=100;y0=50;width=300*(360/96);height=300; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/96) 3 ]);
print('-dpng', '-r400', ['z300_ln_nov_51_20']);


close all;figure();h=pcolor(lon1,lat1,(squeeze(mean(s_dt_en(:,:,11,:),4))-squeeze(mean(s_dt_ln(:,:,11,:),4)))/(100));shading flat;
colorbar;colormap(jet(40));caxis([-3 3]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);

for i=1:720
    for j=1:211
        [r,p] = corrcoef(pr1_nd_enso(:,1), squeeze(z_dt_enso(i,j,1,11,:)));
        if p(1,2)<=0.1;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
end
corr_sst11=corr_sst(1:361,:);corr_sst12=corr_sst(362:720,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
figure();h=pcolor(lon14,lat14,corr_sst2);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,corr_sst2,[-0.4 0.4],'color','k','LineWidth',1);
title(['Correlation HMAPR-Z500: ND, [1951-2020] '],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=200*(360/105);height=200; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
print('-dpng', '-r400', ['corr_prz500_enso_lag0_nd']);







%z5, vimd, qfw, z01 sst11
%z5, vimd, z01



