

%1951-2020
%subtrop
rv_dt2=squeeze(nanmean(u_dt,3));
%u_dt2=squeeze(u_dt(:,:,4,:,:));

t11=squeeze(rv_dt2(:,:,3,:));
corr11=u200_0;
for i=1:720
    for j=1:211
        corr12(i,j)=corr11(i,j)*0;
    end
end
for i=1:72
u0(:,:,i)=t11(:,:,i)+corr12(:,:);
end

%-0.7453 0.91 0.5948

sst11=sst_dt(1:180,:,:,:);sst12=sst_dt(181:360,:,:,:);
sst13=cat(1,sst12,sst11);
sst14=squeeze(sst13(:,:,3,:));
corr11=sst_0;
for i=1:360
    for j=1:180
        if corr11(i,j)>0
        corr12(i,j)=corr11(i,j)*0;
        else corr12(i,j)=nan;
    end
    end
end
for i=1:72
sst0(:,:,i)=sst14(:,:,i)+corr12(:,:);
end


tmp11=t1(1:361,:,:);tmp12=t1(362:720,:,:);
t1_1=cat(1,tmp12,tmp11); 
tmp11=t0(1:361,:,:);tmp12=t0(362:720,:,:);
t0_1=cat(1,tmp12,tmp11); 
%sstp0: 160 260 -15 15: 161:261,  70:111
%ssti0:  60 120 -15 10:  61:121,  80:106
%sstp1: 160 250 -15 15: 161:251,  75:106
%ssti1:  40 100 -15 10:  41:101,  80:106
%sstp2: 160 260 -15 15: 161:281,  70:111
%ssti2:  50 100 -10 10:  51:101,  80:101
%tp0:    10  60 -15 20:  21:121, 111:181              
%tn0:    10  60  30 60:  21:121,  31: 91
%tp1:   -10  60 -25 25: 340:480, 101:201            %lon14
%rv0:    10  60  32 40:  21:121,  71: 87
%u0:     10  60  23 35:  21:121,  81:105

sstp_0 =squeeze(nanmean(nanmean(sst0(161:261,75:106,2:71),1),2));% lag 0          
ssti_0 =squeeze(nanmean(nanmean(sst0(61:121,80:106,2:71),1),2));% lag 0
tp_0 =squeeze(nanmean(nanmean(t0(21:121, 111:181,2:71),1),2));% lag 0          
tn_0 =squeeze(nanmean(nanmean(t0(21:121,  31: 91,2:71),1),2));% lag 0
u_0 =squeeze(nanmean(nanmean(u0(21:121,  81:105,2:71),1),2));% lag 0          
rv_0 =squeeze(nanmean(nanmean(rv0(21:121,71: 87,2:71),1),2));% lag 0
sstp_1 =squeeze(nanmean(nanmean(sst1(161:251,  75:106,2:71),1),2));% lag 1          
ssti_1 =squeeze(nanmean(nanmean(sst1(41:101,  80:106,2:71),1),2));% lag 1
tp_1 =squeeze(nanmean(nanmean(t0_1(340:480, 101:201,2:71),1),2));% lag 1          
sstp_2 =squeeze(nanmean(nanmean(sst2(161:261,  75:106,2:71),1),2));% lag 2          
ssti_2 =squeeze(nanmean(nanmean(sst2(51:101,  80:101,2:71),1),2));% lag 2

%PR: {RV0, U0, SSTp0, Tn0, SSTi0, Tp0, SSTi1, SSTp1, SSTp2, Tp1, SSTi2}

%corrcoef
[r1,p1]=corrcoef(pr_0,sstp_0);
rsstp0=r1(1,2);psstp0=p1(1,2);
[r1,p1]=corrcoef(pr_0,sstp_1);
rsstp1=r1(1,2);psstp1=p1(1,2);
[r1,p1]=corrcoef(pr_0,sstp_2);
rsstp2=r1(1,2);psstp2=p1(1,2);
[r1,p1]=corrcoef(pr_0,ssti_0);
rssti0=r1(1,2);pssti0=p1(1,2);
[r1,p1]=corrcoef(pr_0,ssti_1);
rssti1=r1(1,2);pssti1=p1(1,2);
[r1,p1]=corrcoef(pr_0,ssti_2);
rssti2=r1(1,2);pssti2=p1(1,2);
[r1,p1]=corrcoef(pr_0,tp_0);
rtp0=r1(1,2);ptp0=p1(1,2);
[r1,p1]=corrcoef(pr_0,tp_1);
rtp1=r1(1,2);ptp1=p1(1,2);
[r1,p1]=corrcoef(pr_0,tn_0);
rtn0=r1(1,2);ptn0=p1(1,2);
[r1,p1]=corrcoef(pr_0,u_0);
ru0=r1(1,2);pu0=p1(1,2);
[r1,p1]=corrcoef(pr_0,rv_0);
rrv0=r1(1,2);prv0=p1(1,2);


%
%z05 vimd z8500, z03 z01 sst01 sst11

%pr z05 vimd z01

%z05: vimd z8500, z03 z01 sst01 sst11
%z05  z03 sst11
%z01 sst11
%vimd z850

%
for i=1:720
    for j=1:360
        if z_regrid(i,j)<1000;
            pr5(i,j)=nan;
        end
    end
end


%%

%cross validation

%november
pr1     =npr_0([1:50,61:70],1);
sstp1   =nsstp_1([1:50,61:70],1);
vimdi1  =nvimdi_0([1:50,61:70],1);
z1      =nz1_0([1:50,61:70],1);
z5      =nz5_0([1:50,61:70],1);

pr1     =npr_0(1:60,1);
sstp1   =nsstp_1(1:60,1);
vimdi1  =nvimdi_0(1:60,1);
z1      =nz1_0(1:60,1);
z5      =nz5_0(1:60,1);

x=[ sstp1,z1,z5];
mdl_x=fitlm(x,pr1,'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate
x=[ vimdi1];
mdl_x=fitlm(x,pr1,'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate
x=[z5,vimdi1];
mdl_x=fitlm(x,pr1,'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate

%               sst1    z1        z5        vimdi1
%1961-2020              for 1951-1960      
%single         0.021   0.387  -0.634      -0.598   R2 0.53 0.37
%combined                      -0.483      -0.445   R2 0.58

%1951-1960 & 1971-2020  for 1961-1970      
%single         0.087   0.426   -0.622     -0.654   R2 0.54 0.35
%combined                       -0.508     -0.486   R2 0.59

%1951-1970 & 1981-2020  for 1971-1980      
%single         0.057   0.285   -0.651     -0.569   R2 0.53 0.33
%combined                       -0.563     -0.401   R2 0.59

%1951-1980 & 1991-2020  for 1981-1990      
%single         0.280   0.668   -0.509     -0.521   R2 0.53 0.29
%combined                       -0.444     -0.423   R2 0.51

%1951-1990 & 2001-2020  for 1991-2000      
%single         0.036   0.435   -0.579     -0.605   R2 0.57 0.39
%combined                       -0.461     -0.504   R2 0.59

%1951-2000 & 2011-2020  for 2001-2010      
%single         0.033   0.449   -0.610     -0.612   R2 0.51 0.35
%combined                       -0.480     -0.499   R2 0.57

%1951-2010              for 2010-2020
%single         0.008    0.389   -0.610     -0.611   R2 0.51 0.35
%combined                        -0.483     -0.479   R2 0.55



x=[ sstp1,z1,z5];
mdl_x=fitlm(x,pr1,'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate

% Sum of squared residuals
SSR = sum((prd_pr_nov_ewp - npr_0).^2);
% Total sum of squares
TSS = sum(((npr_0 - mean(npr_0)).^2));
% R squared
Rsquared = 1 - SSR/TSS




%march
pr1     =npr_0(1:60,1);
sstp1   =nsstp_0(1:60,1);
rv1     =nrv_0(1:60,1);

pr1     =  npr_0([1:10,21:70],1);
sstp1   =nsstp_0([1:10,21:70],1);
rv1     =  nrv_0([1:10,21:70],1);
x=[sstp1,rv1];
mdl_x=fitlm(x,pr1,'linear','Intercept',false);
mdl_x.Rsquared
mdl_x.Coefficients.Estimate
 
%                       sst1    rv1
%1961-2020              0.463   0.428   for 1951-1960      
%                               R2  0.49        

%                       sst1    rv1
%1951-1960 1971-2020    0.303   0.471     for 1961-1970      
%                               R2  0.42       

%                       sst1    rv1
%1951-1970 1981-2020    0.372   0.404     for 1971-1980      
%                               R2  0.39    

%                       sst1    rv1
%1951-1980 1991-2020    0.334   0.460     for 1981-1990      
%                               R2  0.40    

%                       sst1    rv1
%1951-1990 2001-2020    0.384   0.480     for 1991-2000      
%                               R2  0.46     

%                       sst1    rv1
%1951-2000 2011-2020    0.319   0.443     for 2001-2010      
%                               R2  0.34 

%                       sst1    rv1
%1951-2010              0.379   0.527   for 2011-2020      
%                               R2  0.46        



