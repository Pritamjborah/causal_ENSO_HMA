%clear all;close all;
%sudo mount -t drvfs F: /mnt/f

k21=[1951,1953,1957,1963,1965,1968,1969, 1972,1976,1977,1979,1982,1986,1987, 1991,1994,1997,2002,2004,2006,2009,2014,2015,2018];
k22=k21-1950;
k31=[1954,1955,1964,1970,1971,1973,1974,1975,1983,1984,1988,1995,1998,1999,2000,2007,2010, 2011,2016,2017,2020];
k32=k31-1950;

k11=[1951,1953,1954,1955,1957,1963,1964,1965,1968,1969,1970,1971, 1972,1973,1974,1975,1976,1977,1979,1982,1983,1984,1986,1987,1988, 1991,1994,1995,1997,1998,1999,2000,2002,2004,2006,2007,2009,2010, 2011,2014,2015,2016,2017,2018,2020];
k12=k11-1950;%enln
k41=[1951,1952,1953, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963,1965, 1966 1967 1968,1969, 1972,1976,1977, 1978 1979, 1980, 1981, 1982, 1985, 1986,1987, 1989, 1990, 1991, 1992 1993 1994, 1996, 1997, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009, 2012, 2013, 2014,2015,2018, 2019];
k42=k41-1950; %en+neutal
k51=[1952, 1954,1955, 1956, 1958, 1959, 1960, 1961, 1962, 1964, 1966, 1967, 1970,1971,1973,1974,1975, 1978, 1980, 1981, 1983,1984, 1985, 1988, 1989, 1990, 1992, 1993, 1995, 1996, 1998, 1999, 2000, 2001, 2003, 2005, 2007, 2008, 2010, 2011, 2012, 2013, 2016,2017,2019, 2020];
k52=k51-1950; %ln+neutral

k61=[1952, 1956, 1958, 1959, 1960, 1961, 1962, 1966, 1967, 1978, 1980, 1981, 1985, 1989, 1990, 1992, 1993, 1996, 2001, 2003, 2005, 2008, 2012, 2013,2019];
k62=k61-1950; %neutral

%ENSO Teleconnections towards Winter Precipitation over High Mountain Asia with a Causal Network Approach
%Tropical Pacific Influence on the Winter Precipitation over High Mountain Asia
%ENSO Teleconnections towards Winter Precipitation over High Mountain Asia
k11m=[1953,1955,1956,1958, 1959,1966,1969,1971,1973,1974,1975,1976,1983,1985,1987,1989, 1992,1995,1996,1998,1999,2000,2006,2008,2009,2010, 2011,2012,2015,2016,2018,2019];
k12m=k11m-1950;%enln
k41m=[1951,1952,1953, 1956, 1957, 1958, 1959, 1960, 1961, 1962, 1963,1965, 1966 1967 1968,1969, 1972,1976,1977, 1978 1979, 1980, 1981, 1982, 1985, 1986,1987, 1989, 1990, 1991, 1992 1993 1994, 1996, 1997, 2001, 2002, 2003, 2004, 2005, 2006, 2008, 2009, 2012, 2013, 2014,2015,2018, 2019];
k42m=k41m-1950; %en+neutal
k51m=[1952, 1954,1955, 1956, 1958, 1959, 1960, 1961, 1962, 1964, 1966, 1967, 1970,1971,1973,1974,1975, 1978, 1980, 1981, 1983,1984, 1985, 1988, 1989, 1990, 1992, 1993, 1995, 1996, 1998, 1999, 2000, 2001, 2003, 2005, 2007, 2008, 2010, 2011, 2012, 2013, 2016,2017,2019, 2020];
k52m=k51m-1950; %ln+neutral
for i=1:32
    pr_enln(:,i)=pr2(:,k12m(1,i)+1);
end
for i=1:32
    sst_dt_enln(:,:,:,i)=sst_dt(:,:,:,k12m(1,i)+1);
end
for i=1:46
    u_dt_enn(:,:,:,:,i)=u_200(:,:,:,:,k52m(1,i)+1);
end
for i=1:32
    t_dt_enln(:,:,:,:,i)=t_dt(:,:,:,:,k12m(1,i)+1);
end

for i=1:25
    z_dt_ne(:,:,:,:,i)=z_dt(:,:,:,:,k62(1,i)+1);
end
for i=1:25
    pr_ne(:,i)=pr2(:,k62(1,i)+1);
end

for i=1:24
    z_dt_en(:,:,:,:,i)=z_dt(:,:,:,:,k22(1,i)+1);
end
for i=1:24
    pr_en(:,i)=pr2(:,k22(1,i)+1);
end
for i=1:21
    z_dt_ln(:,:,:,:,i)=z_dt(:,:,:,:,k32(1,i)+1);
end
for i=1:21
    pr_ln(:,i)=pr2(:,k32(1,i)+1);
end

for i=1:45
    z_dt_enln(:,:,:,:,i)=z_dt(:,:,:,:,k12(1,i)+1);
end
for i=1:45
    pr_enln(:,i)=pr2(:,k12(1,i)+1);
end
for i=1:45
    sst_dt_enln(:,:,:,i)=sst_dt(:,:,:,k12(1,i)+1);
end


for i=1:49
    z_dt_enn(:,:,:,:,i)=z_dt(:,:,:,:,k42(i)+1);
end
for i=1:46
    z_dt_lnn(:,:,:,:,i)=z_dt(:,:,:,:,k52(i)+1);
end
for i=1:49
    pr_enn(:,i)=pr2(:,k42(1,i)+1);
end
for i=1:46
    pr_lnn(:,i)=pr2(:,k52(1,i)+1);
end

for i=1:49
    sst_dt_enn(:,:,:,i)=sst_dt(:,:,:,k42(i)+1);
end
for i=1:46
    sst_dt_lnn(:,:,:,i)=sst_dt(:,:,:,k52(i)+1);
end


for i=1:45
    vimd_dt_enln(:,:,:,i)=vimd_dt(:,:,:,k12(1,i)+1);
end
for i=1:49
    vimd_dt_enn(:,:,:,i)=vimd_dt(:,:,:,k42(1,i)+1);
end
for i=1:46
    vimd_dt_lnn(:,:,:,i)=vimd_dt(:,:,:,k52(1,i)+1);
end

%for i=1:49
%    pr1_nd_enn(i,1)=pr1_nd(k42(i),1);
%end
%for i=1:46
%    pr1_nd_lnn(i,1)=pr1_nd(k52(i),1);
%end

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
%k42:17 43, k52:14:39
%yr=46;
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
tmp21=squeeze(pr_hma_nd(i,j,:));
tmp22=detrend(tmp21,1);
pr_nd_dt(i,j,:)=tmp22;
end
end

for i=1:41
for j=1:41
tmp21=squeeze(pr_hma_n(i,j,:));
tmp22=detrend(tmp21,1);
pr_n_dt(i,j,:)=tmp22;
end
end


for l=1:70
for i=1:41
for j=1:41
pr_hma_jfm(i,j,l)=pr_hma(i,j,1,l+1)+pr_hma(i,j,2,l+1)+pr_hma(i,j,3,l+1); %(1951-2020)
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

pr21=pr_nd_dt(3:35,25:39,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_nd_n=pr22;
pr21=pr_nd_dt(3:35,12:25,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_nd_s=pr22;

for l=1:12
for i=1:41
for j=1:41
tmp21=squeeze(pr_hma(i,j,l,:));
tmp22=detrend(tmp21,1);
pr_hma_dt(i,j,l,:)=tmp22;
end
end
count=l
end

pr21=pr_d_dt(3:35,25:39,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_n=pr22;
pr21=pr_d_dt(3:35,12:25,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_s=pr22;

pr21=pr_d_dt(3:35,25:39,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_d=pr22;

pr_hma_dt_2=squeeze(pr_hma_dt(3:35, 12:39,:,:));
pr_hma_dt_3=squeeze(nanmean(nanmean(pr_hma_dt_2,1),2));

%yr=50;
for k=1:3
for i=1:360
    for j=1:180
        [r,p] = corrcoef(pr_enln(3,:), squeeze(t_dt_enln_av(i,j,k,:)));
%        [r,p] = corrcoef(pr2(3,1+10*(k-1)+1:yr+10*(k-1)+1), squeeze(sst_dt(i,j,3,1+10*(k-1)+1:yr+10*(k-1)+1)));        
        if p(1,2)<=0.1;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end
%corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 
%corr_sst11=corr_sst(1:180,:,:);corr_sst12=corr_sst(181:360,:,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11);

close all;
for i=1:3
figure();h=pcolor(lons,lats,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lons,lats,corr_sst(:,:,i),[0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR - T, lag-',num2str(3-i),': March'],'FontSize',15);ylim([-30 75]);
%title(['HMAPR-SST lag-0: March ', num2str(1951+(i-1)*10),'-',num2str(1950+yr+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;ylim([-30 75]);
set(gca,'XTick',[-180:30:360]);set(gca,'YTick',[-30:15:90]);
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/110) 2.5 ]);
print('-dpng', '-r400', ['corr_prt_enln_m_',num2str(i)]);
%count=i
end

%yr=50;
%for k=1:3
for i=1:720
    for j=1:211
        [r,p] = corrcoef(pr_enn(3,:), squeeze(u_dt_enn(i,j,:)));
%        [r,p] = corrcoef(pr2(11,1+10*(k-1)+1:yr+10*(k-1)+1), squeeze(vimd_dt(i,j,11,1+10*(k-1)+1:yr+10*(k-1)+1)));
                if p(1,2)<=0.1;
            corr_sst(i,j)=r(1,2);
        else corr_sst(i,j)=nan;
        end
    end
%end
%count=k
end
corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for i=1:3
figure();h=pcolor(lon14,lat14,corr_sst2(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,corr_sst2(:,:,i),[-0.5, -0.4 -0.3 0.3 0.4 0.5],'color','k','LineWidth',1);
title(['HMAPR-RV200, lag-',num2str(3-i),': March'],'FontSize',15);
%title(['HMAPR-QFLux lag-0: November ', num2str(1951+(i-1)*10),'-',num2str(1950+yr+(i-1)*10)],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
set(gca,'XTick',[-180:30:360]);set(gca,'YTick',[-30:15:75]);ylim([-30 75]);
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/110) 3 ]);
print('-dpng', '-r400', ['corr_prrv_enln_m_',num2str(i)]);
%count=i
end


close all;
%for i=1:12
figure();h=pcolor(lon1,lat1,corr_sst);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst,[-0.4, -0.3 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-U200, La Nina: March'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
%set(gca,'XTick',[-180:30:180]);set(gca,'YTick',[-30:15:75]);
x0=100;y0=50;width=250*(720/221);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(720/221) 3 ]);
print('-dpng', '-r400', ['corr_lnn_u']);
%count=i
%end

for i=1:720
    for j=1:360
        if lon_regrid(i,j)<0;
        lonr1(i,j)=lon_regrid(i,j)+360;
        else lonr1(i,j)=lon_regrid(i,j);
        end
    end
end
lonr2=lonr1(1:360,:);lonr3=lonr1(361:720,:);
lon_regrid2=cat(1,lonr3,lonr2);lat_regrid2=lat_regrid;
zr2=z_regrid(1:360,:);zr3=z_regrid(361:720,:);
z_regrid2=cat(1,zr3,zr2);

for i=1:360
    for j=1:180
        if lons(i,j)<0;
        lon11(i,j)=lons(i,j)+360;
        else lon11(i,j)=lons(i,j);
        end
    end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lons1=cat(1,lon13,lon12);lats1=lats;

%for k=1:12
for i=1:360
    for j=1:180
        [r,p] = corrcoef(pr2(11,2:71), squeeze(sst_dt(i,j,11,2:71)));
        if p(1,2)<=0.05;
%         if r(1,2)>=0.2 ;
            corr_sst(i,j)=r(1,2);
%         else if r(1,2)<=-0.2; 
%            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j)=nan;
         end
%         end
    end
end
%count=k
%end
corr_sst11=corr_sst(1:180,:,:);corr_sst12=corr_sst(181:360,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 
%62.25 78.25 29.75 43.35, 
a1=[62.25:0.25:78.25];
a2=[29.75:0.25:43.25];

close all;
figure();h=pcolor(lons1,lats1,corr_sst2);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid2,lat_regrid2,z_regrid2,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lons1,lats1,corr_sst2,[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
hold on;plot(a1,29.75*ones(65,1),'color','r','linewidth',2)
hold on;plot(a1,43.25*ones(65,1),'color','r','linewidth',2)
hold on;plot(62.25*ones(55,1),a2,'color','r','linewidth',2)
hold on;plot(78.25*ones(55,1),a2,'color','r','linewidth',2)
title(['HMAPR-SST, Lag-0'],'FontSize',20);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
set(gca,'XTick',[-180:60:360]);set(gca,'YTick',[-30:30:75]);ylim([-30 75]);
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/110) 2.5 ]);
print('-dpng', '-r400', ['corr_nov_sst_l0']);
print('-djpeg', '-r400', ['corr_nov_sst_l0']);

close all;
for i=1:3
figure();h=pcolor(lons1,lats1,corr_sst2(:,:,i+9));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lons1,lats1,corr_sst2(:,:,i+9),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
title(['HMAPR-SST, Lag-1: ND'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
set(gca,'XTick',[-180:30:360]);set(gca,'YTick',[-30:15:75]);ylim([-30 75]);
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/110) 3 ]);
%print('-dpng', '-r400', ['corr_prsst_lnn_lag1_nd']);
count=i
end

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



%%

%1951-2020
z111=squeeze(z_dt(:,:,3,11,:));
%z112=squeeze(z_dt(:,:,3,10,:));

corr31=z_0;
%corr41=z_1;
for i=1:720
    for j=1:211
        corr32(i,j)=corr31(i,j)*0;
%        corr42(i,j)=corr41(i,j)*0;
    end
end

for i=1:72
z12(:,:,i)=z111(:,:,i)+corr32(:,:);
%z22(:,:,i)=z112(:,:,i)+corr42(:,:);
end
%z11 enso: 145 220,   4 30:291:441  91:143
%z01 enso: 140 270, -10 20:281:541 111:191
%z02     : 275 334,  20 45:551:669  61:111    
%z03     : 280 335,  50 75:561:671   1: 51  
%z04     :   0  40,  40 75:  1: 81   1: 71  
%z05     :  35  90,  30 60: 71:181  31: 91
z11 =squeeze(nanmean(nanmean(z22(291:441, 91:143,:),1),2)); % lag 1
z01 =squeeze(nanmean(nanmean(z12(281:541,111:191,:),1),2));% lag 0
z02 =squeeze(nanmean(nanmean(z12(551:669, 61:111,:),1),2));% lag 0
z03 =squeeze(nanmean(nanmean(z12(561:671,  1: 51,:),1),2));% lag 0
z04 =squeeze(nanmean(nanmean(z12(  1: 81,  1: 71,:),1),2));% lag 0          
z05 =squeeze(nanmean(nanmean(z12( 71:181, 31: 91,:),1),2));% lag 0

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
sst111=sst_dt(1:180,:,:,:);sst121=sst_dt(181:360,:,:,:);
sst2=cat(1,sst121,sst111);

sst01=squeeze(sst2(:,:,11,:));
sst11=squeeze(sst2(:,:,10,:));
sst21=squeeze(sst2(:,:, 9,:));

corr01=sst_0;
corr11=sst_1;
corr21=sst_2;
for i=1:360
    for j=1:180
        corr02(i,j)=corr01(i,j)*0;
        corr12(i,j)=corr11(i,j)*0;
        corr22(i,j)=corr21(i,j)*0;
    end
end

for i=1:72
sst02(:,:,i)=sst01(:,:,i)+corr02(:,:);
sst12(:,:,i)=sst11(:,:,i)+corr12(:,:);
sst22(:,:,i)=sst21(:,:,i)+corr22(:,:);
end
%sst02   : 170 240, -10 10:171:241, 80:100
%sst12   : 165 235, -10 15:166:236, 75:100
%sst22   : 165 230, -10 15:166:231, 75:100

sst03=squeeze(sst02(171:241, 80:100,:));
sst13=squeeze(sst12(166:236, 75:100,:));
sst23=squeeze(sst22(166:231, 75:100,:));

%lont1=lon14(171:241,70:101);
%latt1=lat14(171:241,70:101);

sst04=squeeze((nanmean(nanmean(sst03,1),2)));% lag 0
sst14=squeeze((nanmean(nanmean(sst13,1),2)));% lag 1
sst24=squeeze((nanmean(nanmean(sst23,1),2)));% lag 2


%corrcoef
[r1,p1]=corrcoef(pr1_nd_n,sst04(2:71,1));
rsst04=r1(1,2);psst04=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,sst14(2:71,1));
rsst14=r1(1,2);psst14=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,sst24(2:71,1));
rsst24=r1(1,2);psst24=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z11(2:71,1));
rz11=r1(1,2);pz11=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
[r1,p1]=corrcoef(pr1_nd_n,z05(2:71,1));
rz05=r1(1,2);pz05=p1(1,2);


[r1,p1]=corrcoef(z05(2:71,1),sst04(2:71,1));
rsst04=r1(1,2);psst04=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),sst14(2:71,1));
rsst14=r1(1,2);psst14=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),sst24(2:71,1));
rsst24=r1(1,2);psst24=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z11(2:71,1));
rz11=r1(1,2);pz11=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),pr1_nd_n);
%rz05=r1(1,2);pz05=p1(1,2);


%partial corr 
%pr-sst|z
[r1,p1]=partialcorr([pr1_nd,sst1(2:71,1)],[z1(2:71,1)])

%%
%1951-2020
z111=squeeze(z_dt(:,:,3,11,:));
%z112=squeeze(z_dt(:,:,3,10,:));
%z6: -30 20 10 100
corr31=corr_sst;
%corr41=z_1;
for i=1:720
    for j=1:211
        corr32(i,j)=corr31(i,j)*0;
%        corr42(i,j)=corr41(i,j)*0;
    end
end

for i=1:72
z12(:,:,i)=z111(:,:,i)+corr32(:,:);
%z22(:,:,i)=z112(:,:,i)+corr42(:,:);
end
%z11 enso: 145 220,   4 30:291:441  91:143
%z01 enso: 140 270, -20 20:281:541 111:191
%z02     : 270 335,  23 45:541:671  61:105    
%z03     : 280 340,  50 82:561:681   1: 51
%z03     : 210 330,  50 75:561:681   1: 51
%z04     :   0  40,  40 75:  1: 81   1: 71  
%z05     :  35  80,  26 54: 71:161  43: 99
%z6      :  10 100, -30 20: 21:201 111:211
%z11 =squeeze(nanmean(nanmean(z22(291:441, 91:143,:),1),2)); % lag 1
z01 =squeeze(nanmean(nanmean(z12(281:541,111:191,:),1),2));% lag 0
z02 =squeeze(nanmean(nanmean(z12(541:671, 61:105,:),1),2));% lag 0
z03 =squeeze(nanmean(nanmean(z12(561:681,  1: 51,:),1),2));% lag 0
z04 =squeeze(nanmean(nanmean(z12(  1: 81,  1: 71,:),1),2));% lag 0          
z05 =squeeze(nanmean(nanmean(z12( 71:161, 43: 99,:),1),2));% lag 0
z06 =squeeze(nanmean(nanmean(z12( 21:201, 111:211,:),1),2));% lag 0
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
sst111=sst_dt(1:180,:,:,:);sst121=sst_dt(181:360,:,:,:);
sst2=cat(1,sst121,sst111);

sst01=squeeze(sst2(:,:,11,:));
sst11=squeeze(sst2(:,:,10,:));
%sst21=squeeze(sst2(:,:, 9,:));

corr01=sst_0;
corr11=sst_1;
%corr21=sst_2;
for i=1:360
    for j=1:180
        corr02(i,j)=corr01(i,j)*0;
        corr12(i,j)=corr11(i,j)*0;
%        corr22(i,j)=corr21(i,j)*0;
    end
end

for i=1:72
sst02(:,:,i)=sst01(:,:,i)+corr02(:,:);
sst12(:,:,i)=sst11(:,:,i)+corr12(:,:);
%sst22(:,:,i)=sst21(:,:,i)+corr22(:,:);
end
%sst02   : 170 280, -10 15:171:281, 75:100
%sst12   : 165 250, -10 20:166:251, 70:100

sst03=squeeze(sst02(171:281, 75:100,:));
sst13=squeeze(sst12(166:251, 70:100,:));

%lont1=lon14(171:241,70:101);
%latt1=lat14(171:241,70:101);

sst04=squeeze((nanmean(nanmean(sst03,1),2)));% lag 0
sst14=squeeze((nanmean(nanmean(sst13,1),2)));% lag 1
%sst24=squeeze((nanmean(nanmean(sst23,1),2)));% lag 2


%corrcoef
[r1,p1]=corrcoef(pr1_nd,sst04(2:71,1));
rsst04=r1(1,2);psst04=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,sst14(2:71,1));
rsst14=r1(1,2);psst14=p1(1,2);
%[r1,p1]=corrcoef(pr1_nd_n,sst24(2:71,1));
%rsst24=r1(1,2);psst24=p1(1,2);
%[r1,p1]=corrcoef(pr1_nd_n,z11(2:71,1));
%rz11=r1(1,2);pz11=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
[r1,p1]=corrcoef(pr1_nd,z05(2:71,1));
rz05=r1(1,2);pz05=p1(1,2);


[r1,p1]=corrcoef(z05(2:71,1),sst04(2:71,1));
rsst04=r1(1,2);psst04=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),sst14(2:71,1));
rsst14=r1(1,2);psst14=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),sst24(2:71,1));
%rsst24=r1(1,2);psst24=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),z11(2:71,1));
%rz11=r1(1,2);pz11=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),pr1_nd_n);
%rz05=r1(1,2);pz05=p1(1,2);

%%
%1951-2020
close all;
figure();h=pcolor(lon1,lat1,cpnet_z_0);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,cpnet_z_0,[-0.3, 0.3],'color','k','LineWidth',1);
title(['South HMAPR - SST: '],'FontSize',15);ylim([-30 75]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
%set(gca,'XTick',[-180:30:360]);set(gca,'YTick',[-30:15:90]);
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/80) 2.5 ]);
%print('-dpng', '-r400', ['corrd_prsst1_',num2str(i+10)]);

z111=squeeze(z_dt(:,:,3,11,:));
%z112=squeeze(z_dt(:,:,3,10,:));

corr31=cpnet_z_0;
%corr41=z_1;
for i=1:720
    for j=1:211
        corr32(i,j)=corr31(i,j)*0;
%        corr42(i,j)=corr41(i,j)*0;
    end
end

for i=1:72
z12(:,:,i)=z111(:,:,i)+corr32(:,:);
%z22(:,:,i)=z112(:,:,i)+corr42(:,:);
end

sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11);

sst11=squeeze(sst2(:,:,11,:));
corr11=cpnet_sst_0;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sst012(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

sst11=squeeze(sst2(:,:,11,:));
corr11=cpnet_sst_1;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sst112(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

%z01 enso: 160 280, -20 20:321:561 111:191
%z02     : 290 337,  17 43:581:675  65:117   
%z03     : 275 345,  50 75:551:691   1: 51
%z03     : 240 345,  50 75:481:681   1: 51
%z04     :   5  47,  38 75: 11: 95   1: 75  
%z05     :  37  87,  27 55: 75:175  41: 97
%sst02   : 170 280, -10 10:171:281, 80:101
%sst12   : 165 270, -10 10:166:251, 80:101
z01 =squeeze(nanmean(nanmean(z12(321:561, 111:191,:),1),2));% lag 0
z02 =squeeze(nanmean(nanmean(z12(581:675,  65:117,:),1),2));% lag 0
z03 =squeeze(nanmean(nanmean(z12(481:691,   1: 51,:),1),2));% lag 0
z04 =squeeze(nanmean(nanmean(z12( 11: 95,   1: 75,:),1),2));% lag 0          
z05 =squeeze(nanmean(nanmean(z12( 75:175,  41: 97,:),1),2));% lag 0
sst01 =squeeze(nanmean(nanmean(sst012(171:281, 80:101,:),1),2));% lag 0
sst11 =squeeze(nanmean(nanmean(sst112(166:251, 80:101,:),1),2));% lag 0
pr_hma_dt_4=pr_hma_dt_3(11,:)';

%corrcoef
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),sst01(2:71,1));
rsst01=r1(1,2);psst01=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),sst11(2:71,1));
rsst11=r1(1,2);psst11=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
[r1,p1]=corrcoef(pr_hma_dt_4(2:71,1),z05(2:71,1));
rz05=r1(1,2);pz05=p1(1,2);


[r1,p1]=corrcoef(z05(2:71,1),sst04(2:71,1));
rsst04=r1(1,2);psst04=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),sst14(2:71,1));
rsst14=r1(1,2);psst14=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),sst24(2:71,1));
%rsst24=r1(1,2);psst24=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),z11(2:71,1));
%rz11=r1(1,2);pz11=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z01(2:71,1));
rz01=r1(1,2);pz01=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z02(2:71,1));
rz02=r1(1,2);pz02=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z03(2:71,1));
rz03=r1(1,2);pz03=p1(1,2);
[r1,p1]=corrcoef(z05(2:71,1),z04(2:71,1));
rz04=r1(1,2);pz04=p1(1,2);
%[r1,p1]=corrcoef(z05(2:71,1),pr1_nd_n);
%rz05=r1(1,2);pz05=p1(1,2);

%%
qfn1=double(ncread('slp_sp_sst.nc','sp'));
qfn2=reshape(qfn1,[720,361,12,84]);
qfn3=qfn2(:,:,:,10+1:10+72);
for l=1:12
for i=1:720
for j=1:361
tmp21=squeeze(qfn3(i,j,l,:));
tmp22=detrend(tmp21,1);
sp_dt(i,j,l,:)=tmp22;
end
end
count=l
end

pr21=pr_hma_dt(3:35,12:39,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_w=pr22;
pr21=pr_hma_dt(3:35,25:39,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_n=pr22;
pr21=pr_hma_dt(3:35,12:25,:,:);
pr22=squeeze(nanmean(nanmean(pr21,1),2));
pr1_s=pr22;

for l=1:12
%for k=1:3
for i=1:720
    for j=1:361
        [r,p] = corrcoef(pr1_w(l,2:71), squeeze(qfe_dt(i,j,l,2:71)));
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
for i=1:3
figure();h=pcolor(lon,lat,corr_sst(:,:,i+9));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon,lat,corr_sst(:,:,i+9),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ND HMAPR - VIMD ', num2str(month1{i+9})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;box on;
set(gca,'XTick',[0:30:360]);set(gca,'YTick',[-30:15:75]);ylim([-15 60]);
x0=100;y0=100;width=250*(360/80);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/80) 2.5 ]);
%print('-dpng', '-r400', ['corr_pr_vimd_lag0_',num2str(i)]);
count=i
end

close all;
for i=1:12
figure();h=pcolor(lon1,lat1,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,i),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
%title(['ENSO-GLobal SST lag: Month ', num2str(1951+(i-1)*10),'-',num2str(1990+(i-1)*10)],'FontSize',15);
title(['ND HMAPR - SLP ', num2str(month1{i})],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;box on;
set(gca,'XTick',[0:30:360]);set(gca,'YTick',[-30:15:75]);ylim([-15 60]);
x0=100;y0=100;width=250*(360/80);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/80) 2.5 ]);
%print('-dpng', '-r400', ['corr_pr_slp_',num2str(i)]);
count=i
end

%%tropical pathway
%
close all;
figure();h=pcolor(lon,lat,cpn_vimd_0);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon,lat,cpn_vimd_0,[-0.4 0.4],'color','k','LineWidth',1);
title(['South HMAPR-VIMD, lag-0'],'FontSize',15);ylim([-30 45]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
%set(gca,'XTick',[-180:30:360]);set(gca,'YTick',[-30:15:90]);
x0=100;y0=100;width=250*(360/80);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.5*(360/80) 2.5 ]);
print('-dpng', '-r400', ['corr_nw_prvimd_lag0']);


%1951-2020
z11=squeeze(zl875_dt(:,:,2,11,:));
corr11=z700_0;
for i=1:720
    for j=1:361
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
z700m_0(:,:,i)=z11(:,:,i)+corr12(:,:);
end

vimd1=squeeze(vimd_dt(:,:,10,:));
corr11=vimd_1;
for i=1:720
    for j=1:361
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
vimdm_1(:,:,i)=vimd1(:,:,i)+corr12(:,:);
end

qfe01=squeeze(qfe_dt(:,:,10,:));
corr11=qfe_1;
for i=1:720
    for j=1:361
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
qfem_1(:,:,i)=qfe01(:,:,i)+corr12(:,:);
end

sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11);

sst11=squeeze(sst2(:,:,11,:));
corr11=sst_0;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sstm_0(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

sst11=squeeze(sst2(:,:,10,:));
corr11=cpn_sst_1;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sst112(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

%latid=find(lat(1,:)==43)
%lonid=find(lon(1,:)==36)
%qfw_0    : 65 130 -12 12  :

%qfe_0    : 36  90  20 43  : 
%qfw_0    : 30 130 -15 11  :
%qfw_1    : 35 110 -11  5  : 
%vimdio_0 : 38  85 -15 10  : 
%vimdio_1 : 45  75 -10 10  : 
%vimdip_1 :160 210 -4   8  : 
%sstio_0  : 45  85 -20 10  : 
%sstp_0   :180 280 -10 10  :
%sstp_1   :170 270 -10 10  :
%z7001_0  : 55 165 -40 -3  :
%z7002_0  : 50 105   3 26  :
%z7001_1  : 40 240 -20 35  :
%z8501_0  : 70 180 -40 -5  :
%z8502_0  : 70 105   5 27  :
%z8501_1  : 60 200 -25 15  :
ts_qfe_0  =squeeze(nanmean(nanmean(qfem_0(find(lon(:,1) == 36):find(lon(:,1)== 90),find(lat(1,:)==43):find(lat(1,:)== 20),2:71),1),2));
ts_qfw_0  =squeeze(nanmean(nanmean(qfem_0(find(lon(:,1) == 65):find(lon(:,1)==130),find(lat(1,:)==12):find(lat(1,:)==-12),2:71),1),2));
ts_qfw_1  =squeeze(nanmean(nanmean(qfem_1(find(lon(:,1) == 35):find(lon(:,1)==110),find(lat(1,:)== 5):find(lat(1,:)==-11),2:71),1),2));
ts_vimdi_0=squeeze(nanmean(nanmean(vimdm_0(find(lon(:,1)== 38):find(lon(:,1)== 85),find(lat(1,:)==10):find(lat(1,:)==-15),2:71),1),2));
ts_vimdi_1=squeeze(nanmean(nanmean(vimdm_1(find(lon(:,1)== 45):find(lon(:,1)== 75),find(lat(1,:)==10):find(lat(1,:)==-10),2:71),1),2));
ts_vimdp_1=squeeze(nanmean(nanmean(vimdm_1(find(lon(:,1)==160):find(lon(:,1)==210),find(lat(1,:)== 8):find(lat(1,:)== -4),2:71),1),2));
ts_ssti_0 =squeeze(nanmean(nanmean(sstm_0(find(lons1(:,1) == 45.5):find(lons1(:,1)== 85.5),find(lats1(1,:)==10.5):find(lats1(1,:)==-20.5),2:71),1),2));
ts_sstp_0 =squeeze(nanmean(nanmean(sstm_0(find(lons1(:,1) ==179.5):find(lons1(:,1)==279.5),find(lats1(1,:)==10.5):find(lats1(1,:)==-10.5),2:71),1),2));
ts_sstp_1 =squeeze(nanmean(nanmean(sstm_1(find(lons1(:,1) ==169.5):find(lons1(:,1)==269.5),find(lats1(1,:)==10.5):find(lats1(1,:)==-10.5),2:71),1),2));
ts_z700s_0=squeeze(nanmean(nanmean(z700m_0(find(lon(:,1)== 55):find(lon(:,1)==165),find(lat(1,:)==-3):find(lat(1,:)==-40),2:71),1),2));
ts_z700n_0=squeeze(nanmean(nanmean(z700m_0(find(lon(:,1)== 50):find(lon(:,1)==105),find(lat(1,:)==26):find(lat(1,:)==  3),2:71),1),2));
ts_z700_1 =squeeze(nanmean(nanmean(z700m_1(find(lon(:,1)== 40):find(lon(:,1)==240),find(lat(1,:)==35):find(lat(1,:)==-20),2:71),1),2));
%ts_z850s_0=squeeze(nanmean(nanmean(z850m_0(find(lon(:,1)== 70):find(lon(:,1)==180),find(lat(1,:)==-5):find(lat(1,:)==-40),2:71),1),2));
ts_z850s_0=squeeze(nanmean(nanmean(z850m_0(find(lon(:,1)== 75):find(lon(:,1)==190),find(lat(1,:)==20):find(lat(1,:)==-40),2:71),1),2));
ts_z850n_0=squeeze(nanmean(nanmean(z850m_0(find(lon(:,1)== 70):find(lon(:,1)==105),find(lat(1,:)==27):find(lat(1,:)==  5),2:71),1),2));
ts_z850_1 =squeeze(nanmean(nanmean(z850m_1(find(lon(:,1)== 60):find(lon(:,1)==200),find(lat(1,:)==15):find(lat(1,:)==-25),2:71),1),2));




%qfe_0    : 36  90  15 38   : 73:181  105:151
%qfw_0    : 40 120 -20 10   : 81:241  141:201
%vimdio_0 : 38 75 -15 15    : 77:151  151:211   
%sstio_0 : 45 90 -25 5      : 46: 91   85:116
%sstep_0 : 180 280 -10 10   :181:281   80:101
%sstep_1 : 170 270 -10 10   :171:271   80:101
%z500_0  : 65 105 0 22      :131:211  107:151
%z500_0  : 70 100 5 20      :141:201  111:141
sstp_1 =squeeze(nanmean(nanmean(sst112(171:271, 80 :101,:),1),2));% lag 0
ssti_0 =squeeze(nanmean(nanmean(sst012(46:91,   85 :116,:),1),2));% lag 0          
sstp_0 =squeeze(nanmean(nanmean(sst012(181:281, 80 :101,:),1),2));% lag 0
z500_0 =squeeze(nanmean(nanmean(z12(141:201,  111:141,:),1),2));% lag 0
vimd_0 =squeeze(nanmean(nanmean(vimd12(77:151,  151:211,:),1),2));% lag 0
qfe_0  =squeeze(nanmean(nanmean(qfe12(73:181, 105:151,:),1),2));% lag 0
qfw_0  =squeeze(nanmean(nanmean(qfe12(81:241,  141:201,:),1),2));% lag 0

%corrcoef
[r1,p1]=corrcoef(pr_dt_n(2:71,1),sstp_1(2:71,1));
rsstp1=r1(1,2);psstp1=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),sstp_0(2:71,1));
rsstp0=r1(1,2);psstp0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),ssti_0(2:71,1));
rssti0=r1(1,2);pssti0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),z500_0(2:71,1));
rz0=r1(1,2);pz0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),vimd_0(2:71,1));
rvimd0=r1(1,2);pvimd0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),qfe_0(2:71,1));
rqfe0=r1(1,2);pqfe0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),qfw_0(2:71,1));
rqfw0=r1(1,2);pqfw0=p1(1,2);



%z31=z11(401:681, 1:51,:);
%latz31=lat1(401:681, 1:51,:);
%z32=z31.cosd(latz31);
%z33=(squeeze(sum(sum(z32,1),2)))/(squeeze(sum(sum(cosd(latz31),1),2)));
%z33=(squeeze(mean(mean(z32,1),2)));


for i=1:720
    for j=1:360
        if lon_regrid(i,j)<0;
        lon11(i,j)=lon_regrid(i,j)+360;
        else lon11(i,j)=lon_regrid(i,j);
        end
    end
end
lon12=lon11(1:360,:);lon13=lon11(361:720,:);
lon_rg1=cat(1,lon13,lon12);latrg_1=lat_regrid;
zr12=z_regrid(1:360,:);zr13=z_regrid(361:720,:);
z_rg1=cat(1,zr13,zr12);


for i=1:360
    for j=1:180
        if lons(i,j)<0;
        lon11(i,j)=lon1(i,j)+360;
        else lon11(i,j)=lons(i,j);
        end
    end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lon14=cat(1,lon13,lon12);lat14=lats;
sst111=sst_dt(1:180,:,:,:);sst121=sst_dt(181:360,:,:,:);
sst2=cat(1,sst121,sst111);

sst01=squeeze(sst2(:,:,11,:));
sst11=squeeze(sst2(:,:,10,:));
%sst21=squeeze(sst2(:,:, 9,:));

corr01=sst_0;
corr11=sst_1;
%corr21=sst_2;
for i=1:360
    for j=1:180
        corr02(i,j)=corr01(i,j)*0;
        corr12(i,j)=corr11(i,j)*0;
%        corr22(i,j)=corr21(i,j)*0;
    end
end

for i=1:72
sst02(:,:,i)=sst01(:,:,i)+corr02(:,:);
sst12(:,:,i)=sst11(:,:,i)+corr12(:,:);
%sst22(:,:,i)=sst21(:,:,i)+corr22(:,:);
end
%sst02   : 170 280, -10 15:171:281, 75:100
%sst12   : 165 250, -10 20:166:251, 70:100

sst03=squeeze(sst02(171:281, 75:100,:));
sst13=squeeze(sst12(166:251, 70:100,:));

%lont1=lon14(171:241,70:101);
%latt1=lat14(171:241,70:101);

sst04=squeeze((nanmean(nanmean(sst03,1),2)));% lag 0
sst14=squeeze((nanmean(nanmean(sst13,1),2)));% lag 1
%sst24=squeeze((nanmean(nanmean(sst23,1),2)));% lag 2


%%trop pathway whma
z11=squeeze(z_dt(:,:,1,11,:));
corr11=cpn_z500_0;
for i=1:720
    for j=1:211
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
z12(:,:,i)=z11(:,:,i)+corr12(:,:);
end

vimd11=squeeze(vimd_dt(:,:,11,:));
corr11=cpn_vimd_0;
for i=1:720
    for j=1:361
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
vimd12(:,:,i)=vimd11(:,:,i)+corr12(:,:);
end

qfe11=squeeze(qfe_dt(:,:,11,:));
corr11=cpn_qfe_0;
for i=1:720
    for j=1:361
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
qfe12(:,:,i)=qfe11(:,:,i)+corr12(:,:);
end

sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst2=cat(1,sst12,sst11);

sst11=squeeze(sst2(:,:,11,:));
corr11=cpn_sst_0;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sst012(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

sst11=squeeze(sst2(:,:,10,:));
corr11=cpn_sst_1;
for i=1:360
    for j=1:180
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
sst112(:,:,i)=sst11(:,:,i)+corr12(:,:);
end

%qfe_0    : 36  90  15 38   : 73:181  105:151
%qfw_0    : 40 120 -20 10   : 81:241  141:201
%vimdio_0 : 38 75 -15 15    : 77:151  151:211   
%sstio_0 : 45 90 -25 5      : 46: 91   85:116
%sstep_0 : 180 280 -10 10   :181:281   80:101
%sstep_1 : 170 270 -10 10   :171:271   80:101
%z500_0  : 65 105 0 22      :131:211  107:151
%z500_0  : 70 100 5 20      :141:201  111:141
sstp_1 =squeeze(nanmean(nanmean(sst112(171:271, 80 :101,:),1),2));% lag 0
ssti_0 =squeeze(nanmean(nanmean(sst012(46:91,   85 :116,:),1),2));% lag 0          
sstp_0 =squeeze(nanmean(nanmean(sst012(181:281, 80 :101,:),1),2));% lag 0
z500_0 =squeeze(nanmean(nanmean(z12(141:201,  111:141,:),1),2));% lag 0
vimd_0 =squeeze(nanmean(nanmean(vimd12(77:151,  151:211,:),1),2));% lag 0
qfe_0  =squeeze(nanmean(nanmean(qfe12(73:181, 105:151,:),1),2));% lag 0
qfw_0  =squeeze(nanmean(nanmean(qfe12(81:241,  141:201,:),1),2));% lag 0

%corrcoef
[r1,p1]=corrcoef(pr_dt_n(2:71,1),sstp_1(2:71,1));
rsstp1=r1(1,2);psstp1=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),sstp_0(2:71,1));
rsstp0=r1(1,2);psstp0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),ssti_0(2:71,1));
rssti0=r1(1,2);pssti0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),z500_0(2:71,1));
rz0=r1(1,2);pz0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),vimd_0(2:71,1));
rvimd0=r1(1,2);pvimd0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),qfe_0(2:71,1));
rqfe0=r1(1,2);pqfe0=p1(1,2);
[r1,p1]=corrcoef(pr_dt_n(2:71,1),qfw_0(2:71,1));
rqfw0=r1(1,2);pqfw0=p1(1,2);


%%

for k=1:2
for i=1:33
    for j=1:28
        sd1(i,j,k)=std(pr1(i,j,k,:),0,4);
    end
end
end

for k=1:70
    pr11(:,:,:,k)=pr1(:,:,:,k)./sd1(:,:,:);
end

pr2=squeeze(cat(3,squeeze(pr11(:,:,1,:)),squeeze(pr11(:,:,2,:))));


prn1=squeeze(pr1(:,:,1,:));
prd1=squeeze(pr1(:,:,2,:));


%%suptropical pathway


for i=1:360
for j=1:180
if lons(i,j)<0;
lon11(i,j)=lons(i,j)+360;
else lon11(i,j)=lons(i,j);
end
end
end
lon12=lon11(1:180,:);lon13=lon11(181:360,:);
lons1=cat(1,lon13,lon12);lats1=lats;

for k=1:12
for i=1:360
    for j=1:180
        [r,p] = corrcoef(pr2(10,2:71), squeeze(sst_dt(i,j,k,2:71)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end
%corr_sst11=corr_sst(1:361,:);corr_sst12=corr_sst(362:720,:);
%corr_sst2=cat(1,corr_sst12,corr_sst11); 
corr_sst11=corr_sst(1:180,:,:);corr_sst12=corr_sst(181:360,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for i=1:11
figure();h=pcolor(lons1,lats1,corr_sst2(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lons1,lats1,corr_sst2(:,:,i),[-0.5 -0.3, 0.3 0.5],'color','k','LineWidth',1);
title(['HMAPR-SST, Lag-',num2str(11-i),' : October'],'FontSize',15);ylim([-30 60]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/95);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/95) 3 ]);
%print('-dpng', '-r400', ['corr_prsst_oct_l',num2str(10-i)]);
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

rv_dt2=squeeze(nanmean(rv_dt,3));

for k=1:12
for i=1:720
    for j=1:211
        [r,p] = corrcoef(pr2(3,2:71), squeeze(rv_dt2(i,j,3,2:71)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end
corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for k=1:12
figure();h=pcolor(lon1,lat1,corr_sst(:,:,k));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,k),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
title(['HMAPR-OLR, Lag-',num2str(14-k),' : February'],'FontSize',15);ylim([-30 75]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/110) 3 ]);
%print('-dpng', '-r400', ['corr_prolr_feb_l',num2str(14-k)]);
end


for l=1:72
    pass1(:,l)=bandpass(olrjf4(:,l),[1/60, 1/20],1);
end
%lon 481:541 lat 141:161
olrjf3=squeeze(olrjf_dt(481:541,141:161,:,:));

pr_hma_jf=squeeze(nanmean(pr_hma_dt(:,:,1:2,:),3));

for k=1:31
for i=1:41
    for j=1:41
        [r,p] = corrcoef(nanmean(olrjf4(k+1:k+29,2:71),1), pr_hma_jf(i,j,2:71));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end

for k=1:31
for i=1:720
    for j=1:211
        [r,p] = corrcoef(pr2(2,2:71), squeeze(nanmean(olrjf_dt(i,j,k+1:k+29,2:71),3)));
        if p(1,2)<=0.05;
            corr_sst(i,j,k)=r(1,2);
        else corr_sst(i,j,k)=nan;
        end
    end
end
count=k
end

corr_sst11=corr_sst(1:361,:,:);corr_sst12=corr_sst(362:720,:,:);
corr_sst2=cat(1,corr_sst12,corr_sst11); 

close all;
for i=1:31
figure();h=pcolor(lon_hma,lat_hma,corr_sst(:,:,i));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon_hma,lat_hma,corr_sst(:,:,i),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
title(['HMAPR-OLR, Lag-',num2str((31-i)),': February'],'FontSize',15);%ylim([-30 75]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=350*(6/5);height=350; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/110) 3 ]);
%print('-dpng', '-r400', ['corr_prolrd_feb_l',num2str((31-i)/7)]);
end


close all;
for k=1:3
figure();h=pcolor(lon14,lat14,corr_sst2(:,:,k));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
%hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon14,lat14,corr_sst2(:,:,k),[-0.5 -0.3 0.3 0.5],'color','k','LineWidth',1);
title(['HMAPR-U200, Lag-',num2str(3-k),' : March'],'FontSize',15);ylim([-30 75]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/110);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/110) 3 ]);
%print('-dpng', '-r400', ['corr_pru_mar_l',num2str(3-k)]);
end

close all;
for k=1:10
figure();h=pcolor(lon1,lat1,corr_sst(:,:,k));shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,corr_sst(:,:,k),[-0.4 -0.3, 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-VIMD: ',month1{k}],'FontSize',15);ylim([-30 45]);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/80);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/80) 3 ]);
print('-dpng', '-r400', ['corr_prvimd_nov_lag1']);
end



close all;
figure();h=pcolor(lons1,lats1,sst_0);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lons1,lats1,sst_0,[-0.4 -0.3 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-T: December'],'FontSize',15);ylim([-30 75]);;
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_prtav_dec_lag0']);

close all;
figure();h=pcolor(lon1,lat1,u200_0);shading flat;colorbar;colormap(jet(40));caxis([-1 1]);
hold on;load coastlines.mat;hold on;plot(longitude,latitude,'k','linewidth',1);
%hold on;load coastlines.mat;hold on;plot(coastlon,coastlat,'k','linewidth',1);
hold on;[c,h]=contour(lon_regrid,lat_regrid,z_regrid,[1000, 1000],'color',[0.6,0.6,0.6],'LineWidth',1);
hold on;[c,h]=contour(lon1,lat1,u200_0,[-0.4 -0.3 0.3 0.4],'color','k','LineWidth',1);
title(['HMAPR-T: December'],'FontSize',15);
set(gca,'fontsize',15);set(0,'defaultfigurecolor',[1,1,1]);grid on;
x0=100;y0=50;width=250*(360/105);height=250; box on;
set(gcf,'units','points','position',[x0,y0,width,height]);%daspect([4 4 4]);
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 3*(360/105) 3 ]);
%print('-dpng', '-r400', ['corr_prtav_dec_lag0']);


prm1=nanmean(pr_hma_dt(:,:,:,2:71),4);
prsd1=std(pr_hma_dt(:,:,:,2:71),0,4);

prcv1=prsd1/prm1;


%%
x=[ ts_sstp_0/std(ts_sstp_0), ts_rv_0/std(ts_rv_0) ];
mdl_x=fitlm(x,ts_pr_0/std(ts_pr_0),'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate

sst01n=sst01/std(sst01);
sst11n=sst11/std(sst11);
z01n  =z01/std(z01);
z02n  =z02/std(z02);
z03n  =z03/std(z03);
z04n  =z04/std(z04);
z05n  =z05/std(z05);
pr0n  =pr_hma_dt_4/std(pr_hma_dt_4);

%pr, z05 z01 sst11: -0.60 0.42 0.02
%z05, z04 z03 : -0.44 0.23
%z04, z02: 0.64
%z03, z02: -0.57
%z02,z01: .57
%z01, sst11: 0.87




