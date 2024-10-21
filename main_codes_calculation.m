%--------------------------------------------------------------------------
% These codes were written to analyze the change in burned areas in
% Eastern Siberia from 2001 to 2021 and to identify the role of the subpolar
% North Atlantic decadal variability in it. All the codes were developed in
% MATLAB language.
%--------------------------------------------------------------------------
clear all
clc

lon = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lon');
lat = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lat');
yr  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','year');
ba  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','burn_area');
ta  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','total_area');

ba  = 100*ba./ta;
ba  = squeeze(sum(ba,3));

ba_2019 = ba(:,:,19);
ba_2021 = ba(:,:,21);

ba_diff = mean(ba(:,:,11:21),3)-mean(ba(:,:,1:11),3); % 2011-2021 vs 2001-2011

for j = 1:size(ba,2)
    for i = 1:size(ba,1)
        x = squeeze(ba(i,j,1:11));
        y = squeeze(ba(i,j,11:21));
        if (sum(x)==0&sum(y)==0)|sum(x-y)==0
        ba_diff_h(i,j) = 0;
        else
        [h,p] = ttest2(x,y,'Alpha',0.05);
        ba_diff_h(i,j) = h;
        end
    end
end

lon2d = repmat(lon,[1 length(lat)]);
lat2d = repmat(lat,[1 length(lon)])';

lon1d = reshape(lon2d,[144*16 1]);
lat1d = reshape(lat2d,[144*16 1]);

ba_diff_h1d = reshape(ba_diff_h,[144*16 1]);

lonsig = lon1d(find(ba_diff_h1d==1));
latsig = lat1d(find(ba_diff_h1d==1));

myncid  = netcdf.create('ext_fig1abc.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 144);
dimid2  = netcdf.defDim(myncid, 'lat', 16);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));

varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'ba_2019', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'ba_2021', 'double', [dimid1 dimid2]);
varid5  = netcdf.defVar(myncid, 'ba_diff', 'double', [dimid1 dimid2]);
varid6  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid7  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);

netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, ba_2019);
netcdf.putVar(myncid, varid4, ba_2021);
netcdf.putVar(myncid, varid5, ba_diff);
netcdf.putVar(myncid, varid6, lonsig);
netcdf.putVar(myncid, varid7, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lon = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','lon');
lat = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','lat');
pre = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','pre');
tmx = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.tmx.dat.nc','tmx');

pre = pre(:,:,1201:1452);
tmx = tmx(:,:,1201:1452);

pre = reshape(pre,[720*360 252]);
tmx = reshape(tmx,[720*360 252]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);
tmx = tmx(lct,:);

pre_avg = mean(pre.*repmat(weights,[1 252]),1);
tmx_avg = mean(tmx.*repmat(weights,[1 252]),1);

pre_avg  = reshape(pre_avg,[12 21]);
tmx_avg  = reshape(tmx_avg,[12 21]);

pre_clim = mean(pre_avg,2);
tmx_clim = mean(tmx_avg,2);

pre_ann = sum(pre_avg(4:9,:),1);
tmx_ann = mean(tmx_avg(4:9,:),1);

clear lat lat1d lat2d lon lon1d lon2d lct pre tmx weight weights

ndvi = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\ndvi_2001_2022_north_earth.nc','ndvi');
lat  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\ndvi_2001_2022_north_earth.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\ndvi_2001_2022_north_earth.nc','lon');

ndvi   = ndvi(:,:,:,1:21);
ndvi   = reshape(ndvi,[144*16 12*21]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[144 1]);
weight = reshape(weight,[144*16 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[144*16 1]);
lat1d  = reshape(lat2d,[144*16 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

ndvi = ndvi(lct,:);
ndvi_avg  = mean(ndvi.*repmat(weights,[1 252]),1);
ndvi_avg  = reshape(ndvi_avg,[12 21]);
ndvi_clim = mean(ndvi_avg,2);
ndvi_ann  = mean(ndvi_avg(4:9,:),1);

lon = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lon');
lat = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lat');
yr  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','year');
ba  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','burn_area');
ta  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','total_area');

ba  = 100*ba./ta;
ba  = ba(:,:,:,1:21);
ba  = reshape(ba,[144*16 12*21]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[144 1]);
weight = reshape(weight,[144*16 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[144*16 1]);
lat1d  = reshape(lat2d,[144*16 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

ba = ba(lct,:);
ba_avg = mean(ba.*repmat(weights,[1 252]),1);
ba_avg = reshape(ba_avg,[12 21]);
ba_clim = mean(ba_avg(:,1:21),2);
ba_ann  = sum(ba_avg(4:9,:),1);

myncid  = netcdf.create('ext_fig1de.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'mon', 12);
dimid2  = netcdf.defDim(myncid, 'yr', 21);
varid1  = netcdf.defVar(myncid, 'ba_clim', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_clim', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'tmx_clim', 'double', [dimid1]);
varid4  = netcdf.defVar(myncid, 'ndvi_clim', 'double', [dimid1]);
varid5  = netcdf.defVar(myncid, 'ba_ann', 'double', [dimid2]);
varid6  = netcdf.defVar(myncid, 'pre_ann', 'double', [dimid2]);
varid7  = netcdf.defVar(myncid, 'tmx_ann', 'double', [dimid2]);
varid8  = netcdf.defVar(myncid, 'ndvi_ann', 'double', [dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, ba_clim);
netcdf.putVar(myncid, varid2, pre_clim);
netcdf.putVar(myncid, varid3, tmx_clim);
netcdf.putVar(myncid, varid4, ndvi_clim);
netcdf.putVar(myncid, varid5, ba_ann);
netcdf.putVar(myncid, varid6, pre_ann);
netcdf.putVar(myncid, varid7, tmx_ann);
netcdf.putVar(myncid, varid8, ndvi_ann);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri   = 'I:\WORKS\29-Siberian_Fire\ext_fig2\Eastern_siberia_2_2.5degree_2001_2021\';

lon    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lon');
lat    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lat');

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[12*5 1]);
lat1d  = reshape(lat2d,[12*5 1]);

weight = cos(pi.*lat1d/180.0);
weight = weight./mean(weight);

ba     = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'burned_fraction');

pre    = ncread(strcat(diri,'Eastern_siberia_pre_2.5degree_2001_2021.nc'),'pre');
tmx    = ncread(strcat(diri,'Eastern_siberia_tmx_2.5degree_2001_2021.nc'),'tmx');
tmn    = ncread(strcat(diri,'Eastern_siberia_tmn_2.5degree_2001_2021.nc'),'tmn');
vpd    = ncread(strcat(diri,'Eastern_siberia_vpd_2.5degree_2001_2021.nc'),'vpd');
sm     = ncread(strcat(diri,'Eastern_siberia_sm_2.5degree_2001_2021.nc'),'sm');
ws     = ncread(strcat(diri,'Eastern_siberia_ws_2.5degree_2001_2021.nc'),'ws');
frs    = ncread(strcat(diri,'Eastern_siberia_frs_2.5degree_2001_2021.nc'),'frs');

ba     = reshape(sum(ba(:,:,4:9,:),3),[60 21]);
y      = reshape(ba,[60*21 1]);

pre     = reshape(sum(pre(:,:,4:9,:),3),[60 21]);
pre     = (pre-repmat(mean(pre,2),[1 21]))./repmat(std(pre,0,2),[1 21]);
x(:,1)  = reshape(pre,[60*21 1]);

tmx     = reshape(mean(tmx(:,:,4:9,:),3),[60 21]);
tmx     = (tmx-repmat(mean(tmx,2),[1 21]))./repmat(std(tmx,0,2),[1 21]);
x(:,2)  = reshape(tmx,[60*21 1]);

tmn     = reshape(mean(tmn(:,:,4:9,:),3),[60 21]);
tmn     = (tmn-repmat(mean(tmn,2),[1 21]))./repmat(std(tmn,0,2),[1 21]);
x(:,3)  = reshape(tmn,[60*21 1]);

vpd     = reshape(mean(vpd(:,:,4:9,:),3),[60 21]);
vpd     = (vpd-repmat(mean(vpd,2),[1 21]))./repmat(std(vpd,0,2),[1 21]);
x(:,4)  = reshape(vpd,[60*21 1]);

sm      = reshape(mean(sm(:,:,4:9,:),3),[60 21]);
sm      = (sm-repmat(mean(sm,2),[1 21]))./repmat(std(sm,0,2),[1 21]);
x(:,5)  = reshape(sm,[60*21 1]);

ws      = reshape(mean(ws(:,:,4:9,:),3),[60 21]);
ws      = (ws-repmat(mean(ws,2),[1 21]))./repmat(std(ws,0,2),[1 21]);
x(:,6)  = reshape(ws,[60*21 1]);

frs     = reshape(sum(frs(:,:,4:9,:),3),[60 21]);
frs     = (frs-repmat(mean(frs,2),[1 21]))./repmat(std(frs,0,2),[1 21]);
x(:,7)  = reshape(frs,[60*21 1]);

trees = 100;
leaf  = 3;
wuc   = 'on';
Importance = 'on';
rng("default")

y  = reshape(y,[60 21]);
x  = reshape(x,[60 21 7]);

for i = 1:50
train = sort(randsample(21,21*0.8));
test  = setdiff([1:21],train);

x_train = reshape(x(:,train,:),[60*length(train) 7]);
y_train = reshape(y(:,train),[60*length(train) 1]);

x_test = reshape(x(:,test,:),[60*length(test) 7]);
y_test = reshape(y(:,test),[60*length(test) 1]);

net = TreeBagger(trees,x_train,y_train,'OOBPredictorImportance',Importance,'Method','regression',...
                 'OOBPrediction',wuc,'MinLeafSize',leaf);

import(:,i) = net.OOBPermutedPredictorDeltaError;

pred_test = predict(net,x_test);
mab(i)    = mean(abs(pred_test-y_test));
rmse(i)   = sqrt(mean((pred_test-y_test).^2));

clear train test x_train y_train x_test y_test net pred_test
disp(i)
end

myncid  = netcdf.create('ext_fig2a.nc', 'NC_NOCLOBBER');
dimid0  = netcdf.defDim(myncid, 'stat', 3);
dimid1  = netcdf.defDim(myncid, 'factor', 7);
varid1  = netcdf.defVar(myncid, 'import', 'double', [dimid0 dimid1]);
varid2  = netcdf.defVar(myncid, 'mab', 'double', [dimid0]);
varid3  = netcdf.defVar(myncid, 'rmse', 'double', [dimid0]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, [prctile(import,5,2)';prctile(import,50,2)';prctile(import,95,2)']);
netcdf.putVar(myncid, varid2, [prctile(mab,5) prctile(mab,50) prctile(mab,95)]);
netcdf.putVar(myncid, varid3, [prctile(rmse,5) prctile(rmse,50) prctile(rmse,95)]);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri   = 'I:\WORKS\29-Siberian_Fire\ext_fig2\Eastern_siberia_2_2.5degree_2001_2021\';

lon    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lon');
lat    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lat');

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[12*5 1]);
lat1d  = reshape(lat2d,[12*5 1]);

weight = cos(pi.*lat1d/180.0);
weight = weight./mean(weight);

ba     = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'burned_fraction');

pre    = ncread(strcat(diri,'Eastern_siberia_pre_2.5degree_2001_2021.nc'),'pre');
tmx    = ncread(strcat(diri,'Eastern_siberia_tmx_2.5degree_2001_2021.nc'),'tmx');
tmn    = ncread(strcat(diri,'Eastern_siberia_tmn_2.5degree_2001_2021.nc'),'tmn');
vpd    = ncread(strcat(diri,'Eastern_siberia_vpd_2.5degree_2001_2021.nc'),'vpd');
sm     = ncread(strcat(diri,'Eastern_siberia_sm_2.5degree_2001_2021.nc'),'sm');
ws     = ncread(strcat(diri,'Eastern_siberia_ws_2.5degree_2001_2021.nc'),'ws');
frs    = ncread(strcat(diri,'Eastern_siberia_frs_2.5degree_2001_2021.nc'),'frs');

ba     = reshape(sum(ba(:,:,4:9,:),3),[60 21]);
y      = reshape(ba,[60*21 1]);

pre     = reshape(sum(pre(:,:,4:9,:),3),[60 21]);
pre     = (pre-repmat(mean(pre,2),[1 21]))./repmat(std(pre,0,2),[1 21]);
x(:,1)  = reshape(pre,[60*21 1]);

tmx     = reshape(mean(tmx(:,:,4:9,:),3),[60 21]);
tmx     = (tmx-repmat(mean(tmx,2),[1 21]))./repmat(std(tmx,0,2),[1 21]);
x(:,2)  = reshape(tmx,[60*21 1]);

tmn     = reshape(mean(tmn(:,:,4:9,:),3),[60 21]);
tmn     = (tmn-repmat(mean(tmn,2),[1 21]))./repmat(std(tmn,0,2),[1 21]);
x(:,3)  = reshape(tmn,[60*21 1]);

vpd     = reshape(mean(vpd(:,:,4:9,:),3),[60 21]);
vpd     = (vpd-repmat(mean(vpd,2),[1 21]))./repmat(std(vpd,0,2),[1 21]);
x(:,4)  = reshape(vpd,[60*21 1]);

sm      = reshape(mean(sm(:,:,4:9,:),3),[60 21]);
sm      = (sm-repmat(mean(sm,2),[1 21]))./repmat(std(sm,0,2),[1 21]);
x(:,5)  = reshape(sm,[60*21 1]);

ws      = reshape(mean(ws(:,:,4:9,:),3),[60 21]);
ws      = (ws-repmat(mean(ws,2),[1 21]))./repmat(std(ws,0,2),[1 21]);
x(:,6)  = reshape(ws,[60*21 1]);

frs     = reshape(mean(frs(:,:,4:9,:),3),[60 21]);
frs     = (frs-repmat(mean(frs,2),[1 21]))./repmat(std(frs,0,2),[1 21]);
x(:,7)  = reshape(frs,[60*21 1]);

b  = polyfit(x(:,1),x(:,5),1);
x(:,5) = x(:,5) - (x(:,1)*b(1)+b(2));

b  = polyfit(x(:,2),x(:,4),1);
x(:,4) = x(:,4) - (x(:,2)*b(1)+b(2));

trees = 100;
leaf  = 3;
wuc   = 'on';
Importance = 'on';
rng("default")

y  = reshape(y,[60 21]);
x  = reshape(x,[60 21 7]);

for i = 1:50
train = sort(randsample(21,21*0.8));
test  = setdiff([1:21],train);

x_train = reshape(x(:,train,:),[60*length(train) 7]);
y_train = reshape(y(:,train),[60*length(train) 1]);

x_test = reshape(x(:,test,:),[60*length(test) 7]);
y_test = reshape(y(:,test),[60*length(test) 1]);

net = TreeBagger(trees,x_train,y_train,'OOBPredictorImportance',Importance,'Method','regression',...
                 'OOBPrediction',wuc,'MinLeafSize',leaf);

import(:,i) = net.OOBPermutedPredictorDeltaError;

pred_test = predict(net,x_test);
mab(i)    = mean(abs(pred_test-y_test));
rmse(i)   = sqrt(mean((pred_test-y_test).^2));

clear train test x_train y_train x_test y_test net pred_test
disp(i)
end

myncid  = netcdf.create('ext_fig2b.nc', 'NC_NOCLOBBER');
dimid0  = netcdf.defDim(myncid, 'stat', 3);
dimid1  = netcdf.defDim(myncid, 'factor', 7);
varid1  = netcdf.defVar(myncid, 'import', 'double', [dimid0 dimid1]);
varid2  = netcdf.defVar(myncid, 'mab', 'double', [dimid0]);
varid3  = netcdf.defVar(myncid, 'rmse', 'double', [dimid0]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, [prctile(import,5,2)';prctile(import,50,2)';prctile(import,95,2)']);
netcdf.putVar(myncid, varid2, [prctile(mab,5) prctile(mab,50) prctile(mab,95)]);
netcdf.putVar(myncid, varid3, [prctile(rmse,5) prctile(rmse,50) prctile(rmse,95)]);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri   = 'I:\WORKS\29-Siberian_Fire\ext_fig2\Eastern_siberia_2_2.5degree_2001_2021\';

lon    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lon');
lat    = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'lat');

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[12*5 1]);
lat1d  = reshape(lat2d,[12*5 1]);

weight = cos(pi.*lat1d/180.0);
weight = weight./mean(weight);

ba     = ncread(strcat(diri,'Eastern_siberia_burned_fraction_2.5degree_2001_2021.nc'),'burned_fraction');

pre    = ncread(strcat(diri,'Eastern_siberia_pre_2.5degree_2001_2021.nc'),'pre');
tmx    = ncread(strcat(diri,'Eastern_siberia_tmx_2.5degree_2001_2021.nc'),'tmx');
tmn    = ncread(strcat(diri,'Eastern_siberia_tmn_2.5degree_2001_2021.nc'),'tmn');
vpd    = ncread(strcat(diri,'Eastern_siberia_vpd_2.5degree_2001_2021.nc'),'vpd');
sm     = ncread(strcat(diri,'Eastern_siberia_sm_2.5degree_2001_2021.nc'),'sm');
ws     = ncread(strcat(diri,'Eastern_siberia_ws_2.5degree_2001_2021.nc'),'ws');
frs    = ncread(strcat(diri,'Eastern_siberia_frs_2.5degree_2001_2021.nc'),'frs');

ba     = reshape(sum(ba(:,:,4:9,:),3),[60 21]);
y      = reshape(ba,[60*21 1]);

pre     = reshape(sum(pre(:,:,4:9,:),3),[60 21]);
x(:,1)  = reshape(pre,[60*21 1]);

tmx     = reshape(mean(tmx(:,:,4:9,:),3),[60 21]);
x(:,2)  = reshape(tmx,[60*21 1]);

tmn     = reshape(mean(tmn(:,:,4:9,:),3),[60 21]);
x(:,3)  = reshape(tmn,[60*21 1]);

vpd     = reshape(mean(vpd(:,:,4:9,:),3),[60 21]);
x(:,4)  = reshape(vpd,[60*21 1]);

sm      = reshape(mean(sm(:,:,4:9,:),3),[60 21]);
x(:,5)  = reshape(sm,[60*21 1]);

ws      = reshape(mean(ws(:,:,4:9,:),3),[60 21]);
x(:,6)  = reshape(ws,[60*21 1]);

frs     = reshape(mean(frs(:,:,4:9,:),3),[60 21]);
x(:,7)  = reshape(frs,[60*21 1]);

myncid  = netcdf.create('ext_fig2cd.nc', 'NC_NOCLOBBER');
dimid0  = netcdf.defDim(myncid, 'grid', 60*21);
dimid1  = netcdf.defDim(myncid, 'factor', 7);
varid1  = netcdf.defVar(myncid, 'x', 'double', [dimid0 dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, x);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'D:\LENS_CMIP6\FGOALS-g3\pr\';
files_hist   = dir(strcat(diri,'hist_LENs*.nc'));
files_ssp585 = dir(strcat(diri,'ssp585*.nc'));

lat  = ncread(strcat(diri,files_hist(1).name),'lat');
lon  = ncread(strcat(diri,files_hist(1).name),'lon');

[X,Y]  = meshgrid(lat',lon');
clear lat lon

mons = [31 28 31 30 31 30 31 31 30 31 30 31];
pr_sib_ann = zeros(21,110);
pr_nor_ann = zeros(21,110);

landsea = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','landsea');
lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
landsea = reshape(landsea,[180*90 1]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[180 1]);
weight = reshape(weight,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

lct_sib = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights_sib = weight(lct_sib)/mean(weight(lct_sib));

lct_nor = find(lat1d>=45&lat1d<=90&landsea==1);
weights_nor = weight(lct_nor)/mean(weight(lct_nor));

clear weight lon2d lat2d lon1d lat1d landsea
[Xi,Yi]  = meshgrid(lat',lon');

for i = 1:length(files_hist)
    pr_hist   = ncread(strcat(diri,files_hist(i).name),'PRECT')*3600*24*1000;
    pr_ssp585 = ncread(strcat(diri,files_ssp585(i).name),'PRECT')*3600*24*1000;
    pr = cat(3,pr_hist,pr_ssp585);
    clear pr_hist pr_ssp585
    pr        = reshape(pr,[180 80 12 250]);
    pr        = pr(:,:,:,152:172); % 2001-2021
    pr        = reshape(pr,[180 80 12*21]);

    for k=1:size(pr,3)
        pr_2deg(:,:,k) = interp2(X,Y,pr(:,:,k),Xi,Yi);
    end
    clear pr

    pr_2deg(end,:,:) = pr_2deg(end-1,:,:);

    pr  = reshape(pr_2deg,[180*90 12*21]);
    clear pr_2deg
    
    pr_sib  = pr(lct_sib,:);
    pr_sib_avg = mean(pr_sib.*repmat(weights_sib,[1 12*21]),1);
    pr_sib_avg = reshape(pr_sib_avg,[12 21]);
    pr_sib_avg = pr_sib_avg.*repmat(mons',[1 21]);
    pr_sib_ann(:,i) = sum(pr_sib_avg(4:9,:),1)';
    clear pr_sib pr_sib_avg

    pr_nor  = pr(lct_nor,:);
    pr_nor_avg = mean(pr_nor.*repmat(weights_nor,[1 12*21]),1,'omitnan');
    pr_nor_avg = reshape(pr_nor_avg,[12 21]);
    pr_nor_avg = pr_nor_avg.*repmat(mons',[1 21]);
    pr_nor_ann(:,i) = sum(pr_nor_avg(4:9,:),1)';
    clear pr_nor pr_sib_nor pr
    
    disp(i)
end

years = [2001:2021]';

for i = 1:length(files_hist)
[b,bint,r,rint,stats] = regress(pr_sib_ann(:,i),[ones(size(years)) years],0.1);
pr_sib_trd(i) = b(2)*10;

[b,bint,r,rint,stats] = regress(pr_nor_ann(:,i),[ones(size(years)) years],0.1);
pr_nor_trd(i) = b(2)*10;
end

lct_min5  = find(pr_sib_trd<=prctile(pr_sib_trd,4.5));
lct_max5  = find(pr_sib_trd>=prctile(pr_sib_trd,95.5));

myncid  = netcdf.create('fig2a-fgoals-g3.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'mme', 110);
dimid2  = netcdf.defDim(myncid, 'stat', 2);
dimid3  = netcdf.defDim(myncid, 'ext', 5);
varid1  = netcdf.defVar(myncid, 'pr_sib_trd', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pr_sib_trd_stat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pr_nor_trd', 'double', [dimid1]);
varid4  = netcdf.defVar(myncid, 'pr_nor_trd_stat', 'double', [dimid2]);
varid5  = netcdf.defVar(myncid, 'pr_sib_trd_min5', 'double', [dimid3]);
varid6  = netcdf.defVar(myncid, 'pr_sib_trd_max5', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, pr_sib_trd);
netcdf.putVar(myncid, varid2, [prctile(pr_sib_trd,5) prctile(pr_sib_trd,95)]);
netcdf.putVar(myncid, varid3, pr_nor_trd);
netcdf.putVar(myncid, varid4, [prctile(pr_nor_trd,5) prctile(pr_nor_trd,95)]);
netcdf.putVar(myncid, varid5, pr_sib_trd(lct_min5));
netcdf.putVar(myncid, varid6, pr_sib_trd(lct_max5));
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

pre = pre(:,:,1201:1452); % 2001-2021
nansum = sum(isnan(pre),3);
nansum = reshape(nansum,[720*360 1]);

pre = reshape(pre,[720*360 21*12]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct_sib = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125&nansum==0);
lct_nor = find(lat1d>=45&lat1d<=90&nansum==0);
weights_sib = weight(lct_sib)/mean(weight(lct_sib));
weights_nor = weight(lct_nor)/mean(weight(lct_nor));

pre_sib = pre(lct_sib,:);
pre_sib_avg = mean(pre_sib.*repmat(weights_sib,[1 21*12]),1);
pre_sib_avg = reshape(pre_sib_avg,[12 21]);
pre_sib_ann = sum(pre_sib_avg(4:9,:),1);

years = [2001:2021]';
[b,bint,r,rint,stats] = regress(pre_sib_ann',[ones(size(years)) years],0.1)
pre_sib_trd = b(2)*10
CIlower_sib = bint(2,1)*10
CIupper_sib = bint(2,2)*10

pre_nor = pre(lct_nor,:);
pre_nor_avg = mean(pre_nor.*repmat(weights_nor,[1 21*12]),1);
pre_nor_avg = reshape(pre_nor_avg,[12 21]);
pre_nor_ann = sum(pre_nor_avg(4:9,:),1);

years = [2001:2021]';
[b,bint,r,rint,stats] = regress(pre_nor_ann',[ones(size(years)) years],0.1)
pre_nor_trd = b(2)*10
CIlower_nor = bint(2,1)*10
CIupper_nor = bint(2,2)*10

myncid  = netcdf.create('fig2a-obs.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'val1', 1);
dimid2  = netcdf.defDim(myncid, 'val2', 2);
varid1  = netcdf.defVar(myncid, 'pre_sib_trd', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_nor_trd', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pre_sib_trd_CI', 'double', [dimid2]);
varid4  = netcdf.defVar(myncid, 'pre_nor_trd_CI', 'double', [dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, pre_sib_trd);
netcdf.putVar(myncid, varid2, pre_nor_trd);
netcdf.putVar(myncid, varid3, [CIlower_sib CIupper_sib]);
netcdf.putVar(myncid, varid4, [CIlower_nor CIupper_nor]);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

pr_sib_trd = ncread('H:\WORKS\29-Siberian_Fire\figure2\fig2a-fgoals-g3.nc','pr_sib_trd');
lct_min5  = find(pr_sib_trd<=prctile(pr_sib_trd,4.5));
lct_max5  = find(pr_sib_trd>=prctile(pr_sib_trd,95.5));

diri = 'D:\LENS_CMIP6\FGOALS-g3\pr\';
files_hist  = dir(strcat(diri,'hist_LENs*.nc'));
files_rcp85 = dir(strcat(diri,'ssp585*.nc'));

lat  = ncread(strcat(diri,files_hist(1).name),'lat');
lon  = ncread(strcat(diri,files_hist(1).name),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

mons = [31 28 31 30 31 30 31 31 30 31 30 31];

lat     = ncread('land_sea_mask.nc','lat');
lon     = ncread('land_sea_mask.nc','lon');

[Xi,Yi]  = meshgrid(lat',lon');

years = [2001:2021]';

for i = 1:length(files_hist)
    pr_hist   = ncread(strcat(diri,files_hist(i).name),'PRECT')*3600*24*1000;
    pr_rcp85  = ncread(strcat(diri,files_rcp85(i).name),'PRECT')*3600*24*1000;
    pr = cat(3,pr_hist,pr_rcp85);
    clear pr_hist pr_rcp85
    pr        = reshape(pr,[180 80 12 250]);
    pr        = pr(:,:,:,152:172); % 2001-2021
    pr        = reshape(pr,[180 80 12*21]);

    for k=1:size(pr,3)
        pr_2deg(:,:,k) = interp2(X,Y,pr(:,:,k),Xi,Yi);
    end
    clear pr

    pr_2deg(end,:,:) = pr_2deg(end-1,:,:);

    pr = reshape(pr_2deg,[180*90 12 21]);
    clear pr_2deg

    pr = pr.*repmat(mons,[180*90 1 21]);

    pr_ann = squeeze(sum(pr(:,4:9,:),2));

    for j = 1:size(pr_ann,1)
      [b,bint,r,rint,stats] = regress(pr_ann(j,:)',[ones(size(years)) years],0.1);
      pr_trd(j,i) = b(2)*10;
    end

    clear pr pr_ann
    disp(i)
end

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

pr_trd_diff = mean(pr_trd(:,lct_min5),2) - mean(pr_trd(:,lct_max5),2);
pr_trd_diff = reshape(pr_trd_diff,[180 90]);

for i = 1:length(lon1d)
        x = squeeze(pr_trd(i,lct_min5));
        y = squeeze(pr_trd(i,lct_max5));
        [h,p] = ttest(y-x,0,'Alpha',0.5);
        pr_trd_diff_h(i) = h;
end

lonsig = lon1d(find(pr_trd_diff_h==1));
latsig = lat1d(find(pr_trd_diff_h==1));

myncid  = netcdf.create('fig2b.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pr_trd_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, pr_trd_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

lon = lon+180.0;
pre2(1:360,:,:)   = pre(361:720,:,:);
pre2(361:720,:,:) = pre(1:360,:,:);
pre = pre2;
clear pre2

pre = pre(:,:,1201:1452); % 2001-2021
pre = reshape(pre,[720*360 21*12]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');

for j=1:length(lat)
for i=1:length(lon)
    lct = find(lat1d>=lat(j)-1&lat1d<lat(j)+1&lon1d>=lon(i)-1&lon1d<lon(i)+1);
    pre_2deg(i,j,:) = mean(pre(lct,:),1);
end
disp(j)
end

pre = reshape(pre_2deg,[180 90 12*21]); % 1960-2021

nansum = sum(isnan(pre),3);

pre = reshape(pre,[180 90 12 21]);
pre_ann = squeeze(sum(pre(:,:,4:9,:),3));

year  = [2001:2021]';

for j = 1:size(pre_ann,2)
for i = 1:size(pre_ann,1)
      if nansum(i,j) ==0
      [b,bint,r,rint,stats] = regress(squeeze(pre_ann(i,j,:)),[ones(size(year)) year]);
      pre_trd(i,j) = b(2)*10;
      pre_trd_pval(i,j) = stats(3);
      else
      pre_trd(i,j) = -999;
      pre_trd_pval(i,j) = 1;
      end
end
end

pre_trd_pval = reshape(pre_trd_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);
lonsig = lon1d(find(pre_trd_pval<=0.05));
latsig = lat1d(find(pre_trd_pval<=0.05));

myncid  = netcdf.create('fig2c.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pre_trd', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, pre_trd);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

pr_sib_trd = ncread('H:\WORKS\29-Siberian_Fire\figure2\fig2a-fgoals-g3.nc','pr_sib_trd');
lct_min5  = find(pr_sib_trd<=prctile(pr_sib_trd,4.5));
lct_max5  = find(pr_sib_trd>=prctile(pr_sib_trd,95.5));

diri = 'D:\LENS_CMIP6\FGOALS-g3\Z500\';
files_hist  = dir(strcat(diri,'hist_*.nc'));
files_ssp585 = dir(strcat(diri,'ssp585*.nc'));

lat  = ncread(strcat(diri,files_hist(1).name),'lat');
lon  = ncread(strcat(diri,files_hist(1).name),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');

[Xi,Yi]  = meshgrid(lat',lon');

years = [2001:2021]';

for i = 1:length(files_hist)
    zg_hist    = ncread(strcat(diri,files_hist(i).name),'Z500');
    zg_ssp585  = ncread(strcat(diri,files_ssp585(i).name),'Z500');
    zg = cat(3,zg_hist,zg_ssp585);
    clear zg_hist zg_ssp585

    zg  = reshape(zg,[180 80 12 250]);
    zg  = zg(:,:,:,152:172); % 2001-2021
    zg  = reshape(zg,[180 80 12*21]);

    for k=1:size(zg,3)
        zg_2deg(:,:,k) = interp2(X,Y,zg(:,:,k),Xi,Yi);
    end
    clear zg

    zg_2deg(end,:,:) = zg_2deg(end-1,:,:);

    zg = reshape(zg_2deg,[180*90 12 21]);
    clear zg_2deg

    zg_ann = squeeze(mean(zg(:,4:9,:),2));

    for j = 1:size(zg_ann,1)
      [b,bint,r,rint,stats] = regress(zg_ann(j,:)',[ones(size(years)) years],0.1);
      zg_trd(j,i) = b(2)*10;
    end

    clear zg zg_ann
    disp(i)
end

zg_trd_diff = mean(zg_trd(:,lct_min5),2) - mean(zg_trd(:,lct_max5),2);
zg_trd_diff = reshape(zg_trd_diff,[180 90]);

for i = 1:length(lon1d)
        x = squeeze(zg_trd(i,lct_min5));
        y = squeeze(zg_trd(i,lct_max5));
        [h,p] = ttest(y-x,0,'Alpha',0.05);
        zg_trd_diff_h(i) = h;
end

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);
lonsig = lon1d(find(zg_trd_diff_h==1));
latsig = lat1d(find(zg_trd_diff_h==1));

myncid  = netcdf.create('fig2d.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', dimid1);
varid2  = netcdf.defVar(myncid, 'lat', 'double', dimid2);
varid3  = netcdf.defVar(myncid, 'zg_trd_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', dimid3);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', dimid3);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, zg_trd_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'H:\WORKS\29-Siberian_Fire\figure4\';
plev   = ncread(strcat(diri,'hgt.mon.mean.nc'),'level');
zg   = ncread(strcat(diri,'hgt.mon.mean.nc'),'hgt');
zg   = squeeze(zg(:,:,6,1+12*53:12*74)); % 2001-2021

lat  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lat');
lon  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for k=1:size(zg,3)
    zg_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),zg(:,end:-1:1,k),Xi,Yi);
    disp(k)
end
clear zg

zg_2deg(end,:,:) = zg_2deg(end-1,:,:);
zg = reshape(zg_2deg,[180 90 12 21]); % 2001-2021
zg = zg - repmat(mean(zg,1),[180 1 1 1]); % eddy geopotential height
zg_ann = squeeze(mean(zg(:,:,4:9,:),3)); % 2001-2021

year  = [2001:2021]';

for j = 1:size(zg_ann,2)
for i = 1:size(zg_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(zg_ann(i,j,:)),[ones(size(year)) year]);
      zg_trd(i,j) = b(2)*10;
      zg_trd_pval(i,j) = stats(3);
end
end

zg_trd_pval = reshape(zg_trd_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);
lonsig = lon1d(find(zg_trd_pval<=0.05));
latsig = lat1d(find(zg_trd_pval<=0.05));

myncid  = netcdf.create('fig2e.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'zg_trd', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, zg_trd);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

pr_sib_trd = ncread('H:\WORKS\29-Siberian_Fire\figure2\fig2a-fgoals-g3.nc','pr_sib_trd');
lct_min5  = find(pr_sib_trd<=prctile(pr_sib_trd,4.5));
lct_max5  = find(pr_sib_trd>=prctile(pr_sib_trd,95.5));

diri = 'D:\LENS_CMIP6\FGOALS-g3\ts\';
files_hist  = dir(strcat(diri,'hist_LENs*.nc'));
files_rcp85 = dir(strcat(diri,'ssp585*.nc'));

lat  = ncread(strcat(diri,files_hist(1).name),'lat');
lon  = ncread(strcat(diri,files_hist(1).name),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
landsea  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','landsea');
landsea  = reshape(landsea,[180*90 1]);

[Xi,Yi]  = meshgrid(lat',lon');

years = [2001:2021]';

for i = 1:length(files_hist)
    ts_hist   = ncread(strcat(diri,files_hist(i).name),'TS');
    ts_rcp85  = ncread(strcat(diri,files_rcp85(i).name),'TS');
    ts = cat(3,ts_hist,ts_rcp85);
    clear ts_hist ts_rcp85
    ts        = reshape(ts,[180 80 12 250]);
    ts        = ts(:,:,:,152:172); % 2001-2021
    ts        = reshape(ts,[180 80 12*21]);

    for k=1:size(ts,3)
        ts_2deg(:,:,k) = interp2(X,Y,ts(:,:,k),Xi,Yi);
    end
    clear ts

    ts_2deg(end,:,:) = ts_2deg(end-1,:,:);

    ts = reshape(ts_2deg,[180*90 12 21]);
    clear ts_2deg

    ts_ann = squeeze(mean(ts(:,4:9,:),2));

    for j = 1:size(ts_ann,1)
      [b,bint,r,rint,stats] = regress(ts_ann(j,:)',[ones(size(years)) years],0.1);
      ts_trd(j,i) = b(2)*10;
    end

    clear ts ts_ann
    disp(i)
end

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

ts_trd_diff = mean(ts_trd(:,lct_min5),2) - mean(ts_trd(:,lct_max5),2);
ts_trd_diff(find(landsea==1)) = -999;
ts_trd_diff = reshape(ts_trd_diff,[180 90]);

for i = 1:length(lon1d)
        x = squeeze(ts_trd(i,lct_min5));
        y = squeeze(ts_trd(i,lct_max5));
        [h,p] = ttest(y-x,0,'Alpha',0.05);
        ts_trd_diff_h(i) = h;
end

lonsig = lon1d(find(ts_trd_diff_h==1&landsea'==0));
latsig = lat1d(find(ts_trd_diff_h==1&landsea'==0));

myncid  = netcdf.create('fig3a.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'ts_trd_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, ts_trd_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\spna_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\pr_sib_fgoals_g3.mat

for i = 1:110
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5   = find(pr_trd<=prctile(pr_trd,4.5));
lct_max5   = find(pr_trd>=prctile(pr_trd,95.5));

spna_trd_fgoals_g3 = spna_trd([lct_min5 lct_max5]);
pr_trd_fgoals_g3  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\spna_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\pr_sib_accessesm.mat

for i = 1:40
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5   = find(pr_trd<=prctile(pr_trd,12.5));
lct_max5   = find(pr_trd>=prctile(pr_trd,87.5));

spna_trd_accessesm = spna_trd([lct_min5 lct_max5]);
pr_trd_accessesm  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\spna_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\pr_sib_canesm5.mat

for i = 1:50
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5   = find(pr_trd<=prctile(pr_trd,10));
lct_max5   = find(pr_trd>=prctile(pr_trd,90));

spna_trd_canesm5 = spna_trd([lct_min5 lct_max5]);
pr_trd_canesm5  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\spna_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\pr_sib_miroc6.mat

for i = 1:50
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5  = find(pr_trd<=prctile(pr_trd,10));
lct_max5  = find(pr_trd>=prctile(pr_trd,90));

spna_trd_miroc6 = spna_trd([lct_min5 lct_max5]);
pr_trd_miroc6  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\spna_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\pr_sib_mpiesm12.mat

for i = 1:44
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5  = find(pr_trd<=prctile(pr_trd,12));
lct_max5  = find(pr_trd>=prctile(pr_trd,88));

spna_trd_mpiesm12 = spna_trd([lct_min5 lct_max5]);
pr_trd_mpiesm12  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\spna_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\pr_sib_canesm2.mat

for i = 1:50
   [b,a] = polyfit([2001:2021],spna(52:72,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5   = find(pr_trd<=prctile(pr_trd,10));
lct_max5   = find(pr_trd>=prctile(pr_trd,90));

spna_trd_canesm2 = spna_trd([lct_min5 lct_max5]);
pr_trd_canesm2  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\spna_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\pr_sib_mpi_ge.mat

for i = 1:100
   [b,a] = polyfit([2001:2021],spna(152:172,i),1); 
   spna_trd(i) = b(1)*10;
   [b,a] = polyfit([2001:2021],pr_ann_sib(52:72,i),1); 
   pr_trd(i) = b(1)*10;
end

lct_min5  = find(pr_trd<=prctile(pr_trd,5));
lct_max5  = find(pr_trd>=prctile(pr_trd,95));

spna_trd_mpi_ge = spna_trd([lct_min5 lct_max5]);
pr_trd_mpi_ge  = pr_trd([lct_min5 lct_max5]);

clear spna spna_trd pr_ann_sib pr_sib_trd pr_trd lct_min5 lct_max5 ts_ann_res_spna ts_ann_spna

plot([spna_trd_fgoals_g3 spna_trd_accessesm spna_trd_canesm5 spna_trd_miroc6 spna_trd_mpiesm12 spna_trd_canesm2 spna_trd_mpi_ge],...
     [pr_trd_fgoals_g3 pr_trd_accessesm pr_trd_canesm5 pr_trd_miroc6 pr_trd_mpiesm12 pr_trd_canesm2 pr_trd_mpi_ge],'o')

robustfit([spna_trd_fgoals_g3 spna_trd_accessesm spna_trd_canesm5 spna_trd_miroc6 spna_trd_mpiesm12 spna_trd_canesm2 spna_trd_mpi_ge],...
     [pr_trd_fgoals_g3 pr_trd_accessesm pr_trd_canesm5 pr_trd_miroc6 pr_trd_mpiesm12 pr_trd_canesm2 pr_trd_mpi_ge])

[R,P] = corrcoef([spna_trd_fgoals_g3 spna_trd_accessesm spna_trd_canesm5 spna_trd_miroc6 spna_trd_mpiesm12 spna_trd_canesm2 spna_trd_mpi_ge],...
                 [pr_trd_fgoals_g3 pr_trd_accessesm pr_trd_canesm5 pr_trd_miroc6 pr_trd_mpiesm12 pr_trd_canesm2 pr_trd_mpi_ge])

myncid  = netcdf.create('fig3b.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'sample', 70);
varid1  = netcdf.defVar(myncid, 'spna_trd', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pr_trd', 'double', [dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, [spna_trd_fgoals_g3 spna_trd_accessesm spna_trd_canesm2 spna_trd_canesm5 spna_trd_miroc6 spna_trd_mpi_ge spna_trd_mpiesm12]);
netcdf.putVar(myncid, varid2, [pr_trd_fgoals_g3 pr_trd_accessesm pr_trd_canesm2 pr_trd_canesm5 pr_trd_miroc6 pr_trd_mpi_ge pr_trd_mpiesm12]);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lat = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lon');
sst = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','sst'); % since 1854
sst = sst(:,:,1765:2016); % 2001-2021
sst(find(sst<-1e+36)) = NaN;

sst = reshape(sst,[180*89 12 21]);
sst_ann = squeeze(mean(sst(:,4:9,:),2));

years  = [2001:2021]';

for i = 1:size(sst_ann,1)
    nansum = sum(isnan(sst_ann(i,:)));

    if nansum == 0
      [b,bint,r,rint,stats] = regress(sst_ann(i,:)',[ones(size(years)) years]);
      sst_trd(i) = b(2)*10;
      sst_trd_pval(i) = stats(3);
    else
      sst_trd(i) = NaN;
      sst_trd_pval(i) = 1;
    end
end

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*89 1]);
lat1d  = reshape(lat2d,[180*89 1]);
weight = cos(pi*lat/180.0);
weight = repmat(weight',[180 1]);
weight = reshape(weight,[180*89 1]);

lct_glb = find(lat1d>=-45&lat1d<=65&isnan(sst_trd)'==0);
weights = weight(lct_glb)/mean(weight(lct_glb));
sst_trd_glb = mean(sst_trd(lct_glb)'.*weights);

sst_trd = reshape(sst_trd,[180 89])-sst_trd_glb;
sst_trd(isnan(sst_trd)==1) = -999;

lonsig = lon1d(find(sst_trd_pval<=0.05));
latsig = lat1d(find(sst_trd_pval<=0.05));

myncid  = netcdf.create('fig3c.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 89);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'sst_trd', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, sst_trd);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all 
clc

sst = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','sst'); % since 1854
sst = sst(:,:,1273:2016); % 1960-2021
sst(find(sst<0)) = NaN;

sst = reshape(sst,[180*89 12 62]);
sst_ann = squeeze(mean(sst(:,4:9,:),2)); % fire season
%--------------------------------------------------------------------------
diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

pre = pre(:,:,709:1452); % 1960-2021
pre = reshape(pre,[720*360 62*12]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);
pre_avg = mean(pre.*repmat(weights,[1 62*12]),1);
pre_avg  = reshape(pre_avg,[12 62]);
pre_ann = sum(pre_avg(4:9,:),1);
%--------------------------------------------------------------------------
nansum = sum(isnan(sst_ann),2);

for i = 1:size(sst_ann,1)
    if nansum(i)==0
      [b,bint,r,rint,stats] = regress(sst_ann(i,:)',[ones(size(pre_ann')) pre_ann']);
      sst_pre_reg(i) = b(2)*1000;
      sst_pre_pval(i) = stats(3);
    else
      sst_pre_reg(i) = -999;
      sst_pre_pval(i) = 1;
    end
end

sst_pre_reg = reshape(sst_pre_reg,[180 89]);

lat = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lon');
lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*89 1]);
lat1d  = reshape(lat2d,[180*89 1]);
lonsig = lon1d(find(sst_pre_pval<=0.05));
latsig = lat1d(find(sst_pre_pval<=0.05));

myncid  = netcdf.create('fig3d.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 89);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'sst_pre_reg', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, sst_pre_reg);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all 
clc

lat = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lon');
sst = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','sst'); % since 1854
sst = sst(:,:,1273:2016); % 1960-2021
sst(find(sst<-1e+36)) = NaN;

sst = reshape(sst,[180*89 12 62]);
sst = sst - repmat(mean(sst,3),[1 1 62]); % anomaly
sst_ann = squeeze(mean(sst(:,4:9,:),2)); % fire season
nansum = sum(isnan(sst_ann),2);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*89 1]);
lat1d  = reshape(lat2d,[180*89 1]);
weight = cos(pi*lat/180.0);
weight = repmat(weight',[180 1]);
weight = reshape(weight,[180*89 1]);

lct_spna = find(lon1d>=310&lon1d<=350&lat1d>=45&lat1d<=65&nansum==0);
weights_spna = weight(lct_spna)/mean(weight(lct_spna));

lct_glb = find(lat1d>=-45&lat1d<=65&nansum==0);
weights_glb = weight(lct_glb)/mean(weight(lct_glb));
sst_ann_glb = mean(sst_ann(lct_glb,:).*repmat(weights_glb,[1 62]),1);

for i = 1:size(sst_ann,1)  
    if nansum(i) == 0
      [b,bint,r,rint,stats] = regress(sst_ann(i,:)',[ones(size(sst_ann_glb')) sst_ann_glb']);
      sst_ann_res(i,:) = sst_ann(i,:)-(b(2)*sst_ann_glb+b(1));
    else
      sst_ann_res(i,:) = NaN;
    end
end

sst_ann_res_spna = mean(sst_ann_res(lct_spna,:).*repmat(weights_spna,[1 62]),1);
sst_ann_spna = mean(sst_ann(lct_spna,:).*repmat(weights_spna,[1 62]),1);
spna  =  movmean(sst_ann_res_spna,9);

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre  = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

pre = pre(:,:,59*12+1:1452); % 1960-2021
pre = reshape(pre,[720*360 12 62]);
pre = pre- repmat(mean(pre,3),[1 1 62]); % anomaly
pre = reshape(pre,[720*360 62*12]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);

pre_avg = mean(pre.*repmat(weights,[1 62*12]),1);
pre_avg = reshape(pre_avg,[12 62]);
pre_ann = sum(pre_avg(4:9,:),1);

[R,P] = corrcoef(spna,pre_ann)

myncid  = netcdf.create('fig3e.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'yr', 62);
varid1  = netcdf.defVar(myncid, 'spna', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_ann', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pre_ann_movavg', 'double', [dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, spna);
netcdf.putVar(myncid, varid2, pre_ann);
netcdf.putVar(myncid, varid3, movmean(pre_ann,9));
netcdf.close(myncid);
%**************************************************************************
clear all
clc

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

diri = 'H:\WORKS\29-Siberian_Fire\figure4\';
plev   = ncread(strcat(diri,'hgt.mon.mean.nc'),'level');
zg   = ncread(strcat(diri,'hgt.mon.mean.nc'),'hgt');
zg   = squeeze(zg(:,:,6,145:12*74)); % 1960-2021

lat  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lat');
lon  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for k=1:size(zg,3)
    zg_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),zg(:,end:-1:1,k),Xi,Yi);
    disp(k)
end
clear zg

zg_2deg(end,:,:) = zg_2deg(end-1,:,:);

zg = reshape(zg_2deg,[180 90 12 62]); % 1960-2021
zg = zg - repmat(mean(zg,1),[180 1 1 1]); % eddy geopotential height

zg_ann = squeeze(mean(zg(:,:,4:9,:),3));% 1960-2021

for j = 1:size(zg_ann,2)
for i = 1:size(zg_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(zg_ann(i,j,: )),[ones(size(spna)) spna]);
      zg_reg_spna(i,j) = b(2);
      zg_reg_spna_pval(i,j) = stats(3);
end
disp(j)
end

zg_reg_spna_pval = reshape(zg_reg_spna_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);
lonsig = lon1d(find(zg_reg_spna_pval<=0.05));
latsig = lat1d(find(zg_reg_spna_pval<=0.05));

myncid  = netcdf.create('fig4a.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'zg_reg_spna', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, zg_reg_spna);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

diri = 'H:\WORKS\29-Siberian_Fire\figure4\';
plev = ncread(strcat(diri,'uwnd.mon.mean.nc'),'level');
ua   = ncread(strcat(diri,'uwnd.mon.mean.nc'),'uwnd');
ua   = squeeze(ua(:,:,4,145:12*74)); % 1960-2021
va   = ncread(strcat(diri,'vwnd.mon.mean.nc'),'vwnd');
va   = squeeze(va(:,:,4,145:12*74)); % 1960-2021
shum   = ncread(strcat(diri,'shum.mon.mean.nc'),'shum');
shum   = squeeze(shum(:,:,4,145:12*74)); % 1960-2021

lat  = ncread(strcat(diri,'uwnd.mon.mean.nc'),'lat');
lon  = ncread(strcat(diri,'uwnd.mon.mean.nc'),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for k=1:size(ua,3)
    ua_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),ua(:,end:-1:1,k),Xi,Yi);
    va_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),va(:,end:-1:1,k),Xi,Yi);
    shum_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),shum(:,end:-1:1,k),Xi,Yi);
    disp(k)
end
clear zg

ua_2deg(end,:,:) = ua_2deg(end-1,:,:);
ua = reshape(ua_2deg,[180 90 12 62]); % 1960-2021
ua_ann = squeeze(mean(ua(:,:,4:9,:),3));% 1960-2021

va_2deg(end,:,:) = va_2deg(end-1,:,:);
va = reshape(va_2deg,[180 90 12 62]); % 1960-2021
va_ann = squeeze(mean(va(:,:,4:9,:),3));% 1960-2021

shum_2deg(end,:,:) = shum_2deg(end-1,:,:);
shum = reshape(shum_2deg,[180 90 12 62]); % 1960-2021
shum_clim = mean(mean(shum(:,:,4:9,:),4),3);% 1960-2021

for j = 1:size(ua_ann,2)
for i = 1:size(ua_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(ua_ann(i,j,: )),[ones(size(spna)) spna]);
      ua_reg_spna(i,j) = b(2);

      [b,bint,r,rint,stats] = regress(squeeze(va_ann(i,j,: )),[ones(size(spna)) spna]);
      va_reg_spna(i,j) = b(2);
end
disp(j)
end

myncid  = netcdf.create('fig4b.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'ua_reg_spna', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'va_reg_spna', 'double', [dimid1 dimid2]);
varid5  = netcdf.defVar(myncid, 'shum_clim', 'double', [dimid1 dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, ua_reg_spna);
netcdf.putVar(myncid, varid4, va_reg_spna);
netcdf.putVar(myncid, varid5, shum_clim);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

lon = lon+180.0;
pre2(1:360,:,:)   = pre(361:720,:,:);
pre2(361:720,:,:) = pre(1:360,:,:);
pre = pre2;
clear pre2

pre = pre(:,:,709:1452); % 1960-2021
pre = reshape(pre,[720*360 12 62]);
pre = pre- repmat(mean(pre,3),[1 1 62]); % anomaly

pre = reshape(pre,[720*360 62*12]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');

for j=1:length(lat)
for i=1:length(lon)
    lct = find(lat1d>=lat(j)-1&lat1d<lat(j)+1&lon1d>=lon(i)-1&lon1d<lon(i)+1);
    pre_2deg(i,j,:) = mean(pre(lct,:),1);
end
disp(j)
end

pre = reshape(pre_2deg,[180 90 12 62]); % 1960-2021

pre_ann = squeeze(sum(pre(:,:,4:9,:),3));% 1960-2021

for j = 1:size(pre_ann,2)
for i = 1:size(pre_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(pre_ann(i,j,: )),[ones(size(spna)) spna]);
      pre_reg_spna(i,j) = b(2);
      pre_reg_spna_pval(i,j) = stats(3);
end
disp(j)
end

pre_reg_spna_pval = reshape(pre_reg_spna_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

lonsig = lon1d(find(pre_reg_spna_pval<=0.05));
latsig = lat1d(find(pre_reg_spna_pval<=0.05));

pre_reg_spna(find(pre_reg_spna==0)) = -999;

myncid  = netcdf.create('fig4c.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pre_reg_spna', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, pre_reg_spna);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'D:\DCPPC_AMV\CNRM-CM6-1\zg\';
files_neg  = dir(strcat(diri,'zg_Amon_CNRM-CM6-1_dcppC-amv-neg*.nc'));
files_pos  = dir(strcat(diri,'zg_Amon_CNRM-CM6-1_dcppC-amv-pos*.nc'));

lat  = ncread(strcat(diri,files_neg(1).name),'lat');
lon  = ncread(strcat(diri,files_neg(1).name),'lon');
[X,Y] = meshgrid(lat',lon');
clear lat lon

lat  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for i = 1:24
    zg_neg = ncread(strcat(diri,files_neg(i).name),'zg');
    zg_neg = reshape(squeeze(zg_neg(:,:,6,:)),[256 128 12 10]);
    zg_neg = zg_neg - repmat(mean(zg_neg,1),[256 1 1 1]); % eddy geopotential height
    zg_neg = reshape(zg_neg,[256 128 12*10]);

    for k=1:size(zg_neg,3)
        zg_neg_2deg(:,:,k) = interp2(X,Y,zg_neg(:,:,k),Xi,Yi);
    end
    clear zg_neg

    zg_neg_2deg(:,1,:) = zg_neg_2deg(:,2,:);
    zg_neg_2deg(:,end,:) = zg_neg_2deg(:,end-1,:);
    zg_neg_2deg(end,:,:) = zg_neg_2deg(end-1,:,:);

    zg_neg = reshape(zg_neg_2deg,[180*90 12 10]);
    clear zg_neg_2deg

    zg_ann_neg = squeeze(mean(zg_neg(:,4:9,:),2));
    zg_neg_2d(:,:,:,i) = reshape(zg_ann_neg,[180 90 10]);

    clear zg_neg zg_ann_neg
    
    zg_pos = ncread(strcat(diri,files_pos(i).name),'zg');
    zg_pos = reshape(squeeze(zg_pos(:,:,6,:)),[256 128 12 10]);
    zg_pos = zg_pos - repmat(mean(zg_pos,1),[256 1 1 1]); % eddy geopotential height
    zg_pos = reshape(zg_pos,[256 128 12*10]);

    for k=1:size(zg_pos,3)
        zg_pos_2deg(:,:,k) = interp2(X,Y,zg_pos(:,:,k),Xi,Yi);
    end
    clear zg_pos

    zg_pos_2deg(:,1,:) = zg_pos_2deg(:,2,:);
    zg_pos_2deg(:,end,:) = zg_pos_2deg(:,end-1,:);
    zg_pos_2deg(end,:,:) = zg_pos_2deg(end-1,:,:);

    zg_pos = reshape(zg_pos_2deg,[180*90 12 10]);
    clear zg_pos_2deg

    zg_ann_pos = squeeze(mean(zg_pos(:,4:9,:),2));
    zg_pos_2d(:,:,:,i) = reshape(zg_ann_pos,[180 90 10]);

    clear zg_pos zg_ann_pos
    disp(i)
end

zg_diff = mean(mean(zg_pos_2d-zg_neg_2d,4),3);

zg_neg_2d = reshape(zg_neg_2d,[180 90 10*24]);
zg_pos_2d = reshape(zg_pos_2d,[180 90 10*24]);

for j = 1:size(zg_neg_2d,2)
    for i = 1:size(zg_neg_2d,1)
        x = squeeze(zg_neg_2d(i,j,:));
        y = squeeze(zg_pos_2d(i,j,:));
        if (sum(x)==0&sum(y)==0)|sum(x-y)==0
        zg_diff_h(i,j) = 0;
        else
        [h,p] = ttest2(x,y,'Alpha',0.05);
        zg_diff_h(i,j) = h;
        end
    end
end

lon2d = repmat(lon,[1 length(lat)]);
lat2d = repmat(lat,[1 length(lon)])';

lon1d = reshape(lon2d,[180*90 1]);
lat1d = reshape(lat2d,[180*90 1]);

zg_diff_h1d = reshape(zg_diff_h,[180*90 1]);

lonsig = lon1d(find(zg_diff_h1d==1));
latsig = lat1d(find(zg_diff_h1d==1));

myncid  = netcdf.create('fig4d.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'zg_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, zg_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri_hus = 'D:\DCPPC_AMV\CNRM-CM6-1\hus\';
files_pos_hus  = dir(strcat(diri_hus,'hus_Amon_CNRM-CM6-1_dcppC-amv-pos*.nc'));

lat  = ncread(strcat(diri_hus,files_pos_hus(1).name),'lat');
lon  = ncread(strcat(diri_hus,files_pos_hus(1).name),'lon');
[X,Y] = meshgrid(lat',lon');
clear lat lon

lat  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for i = 1:24
    hus_pos = ncread(strcat(diri_hus,files_pos_hus(i).name),'hus')*1000;
    hus_pos = reshape(squeeze(hus_pos(:,:,4,:)),[256 128 12 10]);
    hus_pos = reshape(hus_pos,[256 128 12*10]);

    for k=1:size(hus_pos,3)
        hus_pos_2deg(:,:,k) = interp2(X,Y,hus_pos(:,:,k),Xi,Yi);
    end
    clear hus_pos

    hus_pos_2deg(:,1,:) = hus_pos_2deg(:,2,:);
    hus_pos_2deg(:,end,:) = hus_pos_2deg(:,end-1,:);
    hus_pos_2deg(end,:,:) = hus_pos_2deg(end-1,:,:);

    hus_pos = reshape(hus_pos_2deg,[180*90 12 10]);
    clear hus_pos_2deg

    hus_ann_pos = squeeze(mean(hus_pos(:,4:9,:),2));
    hus_pos_2d(:,:,:,i) = reshape(hus_ann_pos,[180 90 10]);

    clear hus_pos hus_ann_pos
    disp(i)
end

hus_clim = mean(mean(hus_pos_2d,4),3);
hus_clim(isnan(hus_clim)==1) = -999;

myncid  = netcdf.create('fig4e_hus.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'hus_clim', 'double', [dimid1 dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, hus_clim);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri_ua = 'D:\DCPPC_AMV\CNRM-CM6-1\ua\';
files_neg_ua  = dir(strcat(diri_ua,'ua_Amon_CNRM-CM6-1_dcppC-amv-neg*.nc'));
files_pos_ua  = dir(strcat(diri_ua,'ua_Amon_CNRM-CM6-1_dcppC-amv-pos*.nc'));

lat  = ncread(strcat(diri_ua,files_neg_ua(1).name),'lat');
lon  = ncread(strcat(diri_ua,files_neg_ua(1).name),'lon');
[X,Y] = meshgrid(lat',lon');
clear lat lon

lat  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for i = 1:22
    ua_neg = ncread(strcat(diri_ua,files_neg_ua(i).name),'ua');
    ua_neg = reshape(squeeze(ua_neg(:,:,4,:)),[256 128 12 10]);
    ua_neg = reshape(ua_neg,[256 128 12*10]);

    for k=1:size(ua_neg,3)
        ua_neg_2deg(:,:,k) = interp2(X,Y,ua_neg(:,:,k),Xi,Yi);
    end
    clear ua_neg

    ua_neg_2deg(:,1,:) = ua_neg_2deg(:,2,:);
    ua_neg_2deg(:,end,:) = ua_neg_2deg(:,end-1,:);
    ua_neg_2deg(end,:,:) = ua_neg_2deg(end-1,:,:);

    ua_neg = reshape(ua_neg_2deg,[180*90 12 10]);
    clear ua_neg_2deg

    ua_ann_neg = squeeze(mean(ua_neg(:,4:9,:),2));
    ua_neg_2d(:,:,:,i) = reshape(ua_ann_neg,[180 90 10]);

    clear ua_neg ua_ann_neg
    
    ua_pos = ncread(strcat(diri_ua,files_pos_ua(i).name),'ua');
    ua_pos = reshape(squeeze(ua_pos(:,:,4,:)),[256 128 12 10]);
    ua_pos = reshape(ua_pos,[256 128 12*10]);

    for k=1:size(ua_pos,3)
        ua_pos_2deg(:,:,k) = interp2(X,Y,ua_pos(:,:,k),Xi,Yi);
    end
    clear ua_pos

    ua_pos_2deg(:,1,:) = ua_pos_2deg(:,2,:);
    ua_pos_2deg(:,end,:) = ua_pos_2deg(:,end-1,:);
    ua_pos_2deg(end,:,:) = ua_pos_2deg(end-1,:,:);

    ua_pos = reshape(ua_pos_2deg,[180*90 12 10]);
    clear ua_pos_2deg

    ua_ann_pos = squeeze(mean(ua_pos(:,4:9,:),2));
    ua_pos_2d(:,:,:,i) = reshape(ua_ann_pos,[180 90 10]);

    clear ua_pos ua_ann_pos
    disp(i)
end

ua_diff = mean(mean(ua_pos_2d-ua_neg_2d,4),3);
ua_diff(isnan(ua_diff)==1) = -999;

ua_neg_2d = reshape(ua_neg_2d,[180 90 10*22]);
ua_pos_2d = reshape(ua_pos_2d,[180 90 10*22]);

for j = 1:size(ua_neg_2d,2)
    for i = 1:size(ua_neg_2d,1)
        x = squeeze(ua_neg_2d(i,j,:));
        y = squeeze(ua_pos_2d(i,j,:));
        if (sum(x)==0&sum(y)==0)|sum(x-y)==0
        ua_diff_h(i,j) = 0;
        else
        [h,p] = ttest(y-x,0,'Alpha',0.05);
        ua_diff_h(i,j) = h;
        end
    end
end

lon2d = repmat(lon,[1 length(lat)]);
lat2d = repmat(lat,[1 length(lon)])';

lon1d = reshape(lon2d,[180*90 1]);
lat1d = reshape(lat2d,[180*90 1]);

ua_diff_h1d = reshape(ua_diff_h,[180*90 1]);

lonsig = lon1d(find(ua_diff_h1d==1));
latsig = lat1d(find(ua_diff_h1d==1));

myncid  = netcdf.create('fig4e_ua.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'ua_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, ua_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri_va = 'D:\DCPPC_AMV\CNRM-CM6-1\va\';
files_neg_va  = dir(strcat(diri_va,'va_Amon_CNRM-CM6-1_dcppC-amv-neg*.nc'));
files_pos_va  = dir(strcat(diri_va,'va_Amon_CNRM-CM6-1_dcppC-amv-pos*.nc'));

lat  = ncread(strcat(diri_va,files_neg_va(1).name),'lat');
lon  = ncread(strcat(diri_va,files_neg_va(1).name),'lon');
[X,Y] = meshgrid(lat',lon');
clear lat lon

lat  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for i = 1:22
    va_neg = ncread(strcat(diri_va,files_neg_va(i).name),'va');
    va_neg = reshape(squeeze(va_neg(:,:,4,:)),[256 128 12 10]);
    va_neg = reshape(va_neg,[256 128 12*10]);

    for k=1:size(va_neg,3)
        va_neg_2deg(:,:,k) = interp2(X,Y,va_neg(:,:,k),Xi,Yi);
    end
    clear va_neg

    va_neg_2deg(:,1,:) = va_neg_2deg(:,2,:);
    va_neg_2deg(:,end,:) = va_neg_2deg(:,end-1,:);
    va_neg_2deg(end,:,:) = va_neg_2deg(end-1,:,:);

    va_neg = reshape(va_neg_2deg,[180*90 12 10]);
    clear va_neg_2deg

    va_ann_neg = squeeze(mean(va_neg(:,4:9,:),2));
    va_neg_2d(:,:,:,i) = reshape(va_ann_neg,[180 90 10]);

    clear va_neg va_ann_neg
    
    va_pos = ncread(strcat(diri_va,files_pos_va(i).name),'va');
    va_pos = reshape(squeeze(va_pos(:,:,4,:)),[256 128 12 10]);
    va_pos = reshape(va_pos,[256 128 12*10]);

    for k=1:size(va_pos,3)
        va_pos_2deg(:,:,k) = interp2(X,Y,va_pos(:,:,k),Xi,Yi);
    end
    clear va_pos

    va_pos_2deg(:,1,:) = va_pos_2deg(:,2,:);
    va_pos_2deg(:,end,:) = va_pos_2deg(:,end-1,:);
    va_pos_2deg(end,:,:) = va_pos_2deg(end-1,:,:);

    va_pos = reshape(va_pos_2deg,[180*90 12 10]);
    clear va_pos_2deg

    va_ann_pos = squeeze(mean(va_pos(:,4:9,:),2));
    va_pos_2d(:,:,:,i) = reshape(va_ann_pos,[180 90 10]);

    clear va_pos va_ann_pos
    disp(i)
end

va_diff = mean(mean(va_pos_2d-va_neg_2d,4),3);
va_diff(isnan(va_diff)==1) = -999;

va_neg_2d = reshape(va_neg_2d,[180 90 10*22]);
va_pos_2d = reshape(va_pos_2d,[180 90 10*22]);

for j = 1:size(va_neg_2d,2)
    for i = 1:size(va_neg_2d,1)
        x = squeeze(va_neg_2d(i,j,:));
        y = squeeze(va_pos_2d(i,j,:));
        if (sum(x)==0&sum(y)==0)|sum(x-y)==0
        va_diff_h(i,j) = 0;
        else
        [h,p] = ttest(y-x,0,'Alpha',0.05);
        va_diff_h(i,j) = h;
        end
    end
end

lon2d = repmat(lon,[1 length(lat)]);
lat2d = repmat(lat,[1 length(lon)])';

lon1d = reshape(lon2d,[180*90 1]);
lat1d = reshape(lat2d,[180*90 1]);

va_diff_h1d = reshape(va_diff_h,[180*90 1]);

lonsig = lon1d(find(va_diff_h1d==1));
latsig = lat1d(find(va_diff_h1d==1));

myncid  = netcdf.create('fig4e_va.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'va_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, va_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'D:\DCPPC_AMV\CNRM-CM6-1\pr\';
files_neg  = dir(strcat(diri,'pr_Amon_CNRM-CM6-1_dcppC-amv-neg*.nc'));
files_pos  = dir(strcat(diri,'pr_Amon_CNRM-CM6-1_dcppC-amv-pos*.nc'));

lat  = ncread(strcat(diri,files_neg(1).name),'lat');
lon  = ncread(strcat(diri,files_neg(1).name),'lon');
[X,Y] = meshgrid(lat',lon');
clear lat lon

mons = [31 28 31 30 31 30 31 31 30 31 30 31];

lat  = ncread('I:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon  = ncread('I:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi] = meshgrid(lat',lon');

for i = 1:24
    pr_neg   = ncread(strcat(diri,files_neg(i).name),'pr')*3600*24;
    pr_neg = reshape(pr_neg,[256*128 12 10]);
    pr_neg = pr_neg.*repmat(mons,[256*128 1 10]);
    pr_neg = reshape(pr_neg,[256 128 12 10]);

    pr_neg = reshape(pr_neg,[256 128 12*10]);

    for k=1:size(pr_neg,3)
        pr_neg_2deg(:,:,k) = interp2(X,Y,pr_neg(:,:,k),Xi,Yi);
    end
    clear pr_neg

    pr_neg_2deg(:,1,:) = pr_neg_2deg(:,2,:);
    pr_neg_2deg(:,end,:) = pr_neg_2deg(:,end-1,:);
    pr_neg_2deg(end,:,:) = pr_neg_2deg(end-1,:,:);

    pr_neg = reshape(pr_neg_2deg,[180*90 12 10]);
    clear pr_neg_2deg

    pr_ann_neg = squeeze(sum(pr_neg(:,4:9,:),2));
    pr_neg_2d(:,:,:,i) = reshape(pr_ann_neg,[180 90 10]);

    clear pr_neg pr_ann_neg
    
    pr_pos   = ncread(strcat(diri,files_pos(i).name),'pr')*3600*24;
    pr_pos = reshape(pr_pos,[256*128 12 10]);
    pr_pos = pr_pos.*repmat(mons,[256*128 1 10]);
    pr_pos = reshape(pr_pos,[256 128 12 10]);

    pr_pos = reshape(pr_pos,[256 128 12*10]);

    for k=1:size(pr_pos,3)
        pr_pos_2deg(:,:,k) = interp2(X,Y,pr_pos(:,:,k),Xi,Yi);
    end
    clear pr_pos

    pr_pos_2deg(:,1,:) = pr_pos_2deg(:,2,:);
    pr_pos_2deg(:,end,:) = pr_pos_2deg(:,end-1,:);
    pr_pos_2deg(end,:,:) = pr_pos_2deg(end-1,:,:);

    pr_pos = reshape(pr_pos_2deg,[180*90 12 10]);
    clear pr_pos_2deg

    pr_ann_pos = squeeze(sum(pr_pos(:,4:9,:),2));
    pr_pos_2d(:,:,:,i) = reshape(pr_ann_pos,[180 90 10]);

    clear pr_pos pr_ann_pos
    disp(i)
end

pr_diff = mean(mean(pr_pos_2d-pr_neg_2d,4),3);

pr_neg_2d = reshape(pr_neg_2d,[180 90 10*24]);
pr_pos_2d = reshape(pr_pos_2d,[180 90 10*24]);

for j = 1:size(pr_neg_2d,2)
    for i = 1:size(pr_neg_2d,1)
        x = squeeze(pr_neg_2d(i,j,:));
        y = squeeze(pr_pos_2d(i,j,:));
        if (sum(x)==0&sum(y)==0)|sum(x-y)==0
        pr_diff_h(i,j) = 0;
        else
        [h,p] = ttest2(x,y,'Alpha',0.05);
        pr_diff_h(i,j) = h;
        end
    end
end

lon2d = repmat(lon,[1 length(lat)]);
lat2d = repmat(lat,[1 length(lon)])';

lon1d = reshape(lon2d,[180*90 1]);
lat1d = reshape(lat2d,[180*90 1]);

pr_diff_h1d = reshape(pr_diff_h,[180*90 1]);

lonsig = lon1d(find(pr_diff_h1d==1));
latsig = lat1d(find(pr_diff_h1d==1));

myncid  = netcdf.create('fig4f.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pr_diff', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, pr_diff);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\spna_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\pr_sib_fgoals_g3.mat

spna_ens = spna(111:172,:);
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 110]);
pr_sib_ens     = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\spna_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\pr_sib_accessesm.mat

spna_ens = cat(2,spna_ens,spna(111:172,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 40]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\spna_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\pr_sib_canesm5.mat

spna_ens = cat(2,spna_ens,spna(111:172,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\spna_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\pr_sib_miroc6.mat

spna_ens = cat(2,spna_ens,spna(111:172,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\spna_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\pr_sib_mpiesm12.mat

spna_ens = cat(2,spna_ens,spna(111:172,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 44]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\spna_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\pr_sib_canesm2.mat

spna_ens = cat(2,spna_ens,spna(11:72,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\spna_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\pr_sib_mpi_ge.mat

spna_ens = cat(2,spna_ens,spna(111:172,:));
pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 100]);
pr_ann_sib = (pr_ann_sib_int(11:72,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[62 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[62 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib);

clear spna pr_ann_sib pr_ann_sib_int ts_ann_res_spna ts_ann_spna

for i = 1:size(spna_ens,2)
[R,P] = corrcoef(spna_ens(:,i),pr_sib_ens(:,i));
corr(i) = R(2,1);
pval(i) = P(2,1);
end 

lct = find(corr>0&pval<0.1);

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');
spna_ens = spna_ens(:,lct);
pr_sib_ens = pr_sib_ens(:,lct);

for i = 1:size(spna_ens,2)
b = polyfit(spna_ens(:,i),pr_sib_ens(:,i),1);
pr_sib_ens_adj(:,i) = pr_sib_ens(:,i)+b(1)*(spna-spna_ens(:,i));
end

for i = 1:size(spna_ens,2)
b = polyfit([2001:2021]',pr_sib_ens(42:62,i),1);
pr_sib_ens_trd(i) = b(1)*10;
b = polyfit([2001:2021]',pr_sib_ens_adj(42:62,i),1);
pr_sib_ens_trd_adj(i) = b(1)*10;
end

mean(pr_sib_ens_trd)
mean(pr_sib_ens_trd_adj)

median(pr_sib_ens_trd_adj)

mean(pr_sib_ens_trd_adj)-mean(pr_sib_ens_trd)

npts = 20;
pts = linspace(min(pr_sib_ens_trd)-0.2,max(pr_sib_ens_trd)+0.2,npts); % points to evaluate the estimator
[pr_sib_ens_trd_f,pr_sib_ens_trd_xi] = ksdensity(pr_sib_ens_trd,pts);
pts = linspace(min(pr_sib_ens_trd_adj)-0.2,max(pr_sib_ens_trd_adj)+0.2,npts); % points to evaluate the estimator
[pr_sib_ens_trd_adj_f,pr_sib_ens_trd_adj_xi] = ksdensity(pr_sib_ens_trd_adj,pts);

s = RandStream('mlfg6331_64'); 
for i = 1:1000
sample = randsample(s,length(lct),length(lct),true);
spna_pr_sigma(i) = mean(pr_sib_ens_trd_adj(sample))-mean(pr_sib_ens_trd(sample));
disp(i)
end

prctile(spna_pr_sigma,5)
prctile(spna_pr_sigma,95)

myncid  = netcdf.create('fig5a.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'sample', length(lct));
dimid2  = netcdf.defDim(myncid, 'stat', 6);
dimid3  = netcdf.defDim(myncid, 'pdf', npts);
varid0  = netcdf.defVar(myncid, 'lct', 'double', [dimid1]);
varid1  = netcdf.defVar(myncid, 'pr_sib_ens_trd', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pr_sib_ens_trd_adj', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pr_sib_ens_trd_stat', 'double', [dimid2]);
varid4  = netcdf.defVar(myncid, 'pr_sib_ens_trd_adj_stat', 'double', [dimid2]);
varid5  = netcdf.defVar(myncid, 'pr_sib_ens_trd_f', 'double', [dimid3]);
varid6  = netcdf.defVar(myncid, 'pr_sib_ens_trd_xi', 'double', [dimid3]);
varid7  = netcdf.defVar(myncid, 'pr_sib_ens_trd_adj_f', 'double', [dimid3]);
varid8  = netcdf.defVar(myncid, 'pr_sib_ens_trd_adj_xi', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid0, lct);
netcdf.putVar(myncid, varid1, pr_sib_ens_trd);
netcdf.putVar(myncid, varid2, pr_sib_ens_trd_adj);
netcdf.putVar(myncid, varid3, [prctile(pr_sib_ens_trd,5) prctile(pr_sib_ens_trd,25) median(pr_sib_ens_trd)... 
                               prctile(pr_sib_ens_trd,75) prctile(pr_sib_ens_trd,95) mean(pr_sib_ens_trd)]);
netcdf.putVar(myncid, varid4, [prctile(pr_sib_ens_trd_adj,5) prctile(pr_sib_ens_trd_adj,25) median(pr_sib_ens_trd_adj)... 
                               prctile(pr_sib_ens_trd_adj,75) prctile(pr_sib_ens_trd_adj,95) mean(pr_sib_ens_trd_adj)]);
netcdf.putVar(myncid, varid5, pr_sib_ens_trd_f);
netcdf.putVar(myncid, varid6, pr_sib_ens_trd_xi);
netcdf.putVar(myncid, varid7, pr_sib_ens_trd_adj_f);
netcdf.putVar(myncid, varid8, pr_sib_ens_trd_adj_xi);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\spna_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\pr_sib_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\tasmax_sib_fgoals_g3.mat

[lcti lctj] = find(tasmax_ann_sib>1000);
tasmax_ann_sib(lcti,lctj) = (tasmax_ann_sib(lcti-1,lctj)+tasmax_ann_sib(lcti+1,lctj))*0.5; 

spna_ens = spna(101:206,:); % 1950-2055

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 110]);
pr_sib_ens     = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]); % 1950-2055
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = pr_sib_forced;

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 110]);
tasmax_sib_ens     = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = tasmax_sib_forced;

clear spna pr_ann_sib pr_ann_sib_int pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\spna_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\pr_sib_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\tasmax_sib_accessesm.mat

spna_ens = cat(2,spna_ens,spna(101:206,:)); % 1950-2055

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 40]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 40]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\spna_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\pr_sib_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\tasmax_sib_canesm5.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\spna_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\pr_sib_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\tasmax_sib_miroc6.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\spna_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\pr_sib_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\tasmax_sib_mpiesm12.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 44]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 44]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\spna_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\pr_sib_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\tasmax_sib_canesm2.mat

spna_ens = cat(2,spna_ens,spna(1:106,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\spna_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\pr_sib_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\tas_sib_mpi_ge.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 100]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tas_ann_sib - repmat(mean(tas_ann_sib,2),[1 100]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tas_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tas_ann_sib(1:106,:)-repmat(mean(tas_ann_sib(52:72,:),1),[106 1]))./repmat(std(tas_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

lct = ncread('H:\WORKS\29-Siberian_Fire\figure5\fig5a.nc','lct');
spna_ens = spna_ens(:,lct);

b0 = polyfit([2001:2021]',spna(42:62),1);
spna_trd = b0(1)*10;

for i = 1:length(lct)
for j = 1:52
    b1 = polyfit([2021:2050]',spna_ens(j+20:j+49,i),1);
    spna_ens_trd(j) = b1(1)*10;
    clear b1
end

lct_max_corr         = find(spna_ens_trd==max(spna_ens_trd))+20;
spna_ens_opt(:,i)     = spna_ens(lct_max_corr-20:lct_max_corr+34,i);
pr_ens_opt(:,i)      = pr_sib_ens(lct_max_corr-20:lct_max_corr+34,i);
tasmax_ens_opt(:,i)  = tasmax_sib_ens(lct_max_corr-20:lct_max_corr+34,i);

clear lct_max_corr
end

pr_sib_ens2 = mean(pr_sib_ens2,2);
tasmax_sib_ens2 = mean(tasmax_sib_ens2,2);

pr_sib_ens2_hist     = pr_sib_ens2(11:72);
tasmax_sib_ens2_hist = tasmax_sib_ens2(11:72);

pre_ann = ncread('obs_pr_tasmax.nc','pre_ann');
tasmax_ann = ncread('obs_pr_tasmax.nc','tmx_ann');

pre_ann = (pre_ann-mean(pre_ann(42:62)))/std(pre_ann(42:62),0,1);
tasmax_ann = (tasmax_ann-mean(tasmax_ann(42:62)))/std(tasmax_ann(42:62),0,1);

scale_pr     = polyfit(pr_sib_ens2_hist,pre_ann,1);
scale_tasmax = polyfit(tasmax_sib_ens2_hist,tasmax_ann,1);

pr_ens_opt = pr_ens_opt + repmat(pr_sib_ens2(52:106)*scale_pr(1),[1 size(pr_ens_opt,2)]); % 2001-2050
tasmax_ens_opt = tasmax_ens_opt + repmat(tasmax_sib_ens2(52:106)*scale_tasmax(1),[1 size(tasmax_ens_opt,2)]);

npts = 20;

for i = 1:7
pr_ens_opt_temp  = reshape(pr_ens_opt(16+(i-1)*5:25+(i-1)*5,:),[1 10*size(pr_ens_opt,2)]);
pr_ens_forced(i) = mean(pr_sib_ens2(67+(i-1)*5:76+(i-1)*5));
pts = linspace(min(pr_ens_opt_temp)-0.2,max(pr_ens_opt_temp)+0.2,npts); % points to evaluate the estimator
[pr_ens_opt_temp_f,pr_ens_opt_temp_xi] = ksdensity(pr_ens_opt_temp,pts);
pr_ens_opt_stat(:,i) = [prctile(pr_ens_opt_temp,5) prctile(pr_ens_opt_temp,25) median(pr_ens_opt_temp)... 
                        prctile(pr_ens_opt_temp,75) prctile(pr_ens_opt_temp,95) mean(pr_ens_opt_temp)]';
pr_ens_opt_f(:,i)    = pr_ens_opt_temp_f';
pr_ens_opt_xi(:,i)   = pr_ens_opt_temp_xi';

clear pr_ens_opt_temp pts pr_ens_opt_temp_f pr_ens_opt_temp_xi

tasmax_ens_opt_temp  = reshape(tasmax_ens_opt(16+(i-1)*5:25+(i-1)*5,:),[1 10*size(tasmax_ens_opt,2)]);
tasmax_ens_forced(i) = mean(tasmax_sib_ens2(67+(i-1)*5:76+(i-1)*5));
pts = linspace(min(tasmax_ens_opt_temp)-0.2,max(tasmax_ens_opt_temp)+0.2,npts); % points to evaluate the estimator
[tasmax_ens_opt_temp_f,tasmax_ens_opt_temp_xi] = ksdensity(tasmax_ens_opt_temp,pts);
tasmax_ens_opt_stat(:,i) = [prctile(tasmax_ens_opt_temp,5) prctile(tasmax_ens_opt_temp,25) median(tasmax_ens_opt_temp)... 
                            prctile(tasmax_ens_opt_temp,75) prctile(tasmax_ens_opt_temp,95) mean(tasmax_ens_opt_temp)]';
tasmax_ens_opt_f(:,i)    = tasmax_ens_opt_temp_f';
tasmax_ens_opt_xi(:,i)   = tasmax_ens_opt_temp_xi';

clear tasmax_ens_opt_temp pts tasmax_ens_opt_temp_f tasmax_ens_opt_temp_xi
end

for i = 1:length(lct)
b_pr     = polyfit([2021:2050],pr_ens_opt(21:50,i),1);
b_tasmax = polyfit([2021:2050],tasmax_ens_opt(21:50,i),1);
ba_frac(i) = b_pr(1)*30*(-0.64)+b_tasmax(1)*30*(0.22);
clear b_pr b_tasmax
end

max(ba_frac)
min(ba_frac)

bins = [-0.6:0.3:1.5];
for i = 1:length(bins)-1
lct2 = find(ba_frac>=bins(i)&ba_frac<=bins(i+1));
ba_frac_prob(i) = length(lct2)/length(lct);
end

myncid  = netcdf.create('fig5bcd.nc', 'NC_NOCLOBBER');
dimid0  = netcdf.defDim(myncid, 'decade', 7);
dimid1  = netcdf.defDim(myncid, 'stat', 6);
dimid2  = netcdf.defDim(myncid, 'npts', npts);
dimid3  = netcdf.defDim(myncid, 'mem', length(lct));
dimid4  = netcdf.defDim(myncid, 'year0', 21);
dimid5  = netcdf.defDim(myncid, 'year1', 50);
dimid6  = netcdf.defDim(myncid, 'bins', length(bins)-1);
varid11  = netcdf.defVar(myncid, 'pr_ens_opt_stat', 'double', [dimid1 dimid0]);
varid12  = netcdf.defVar(myncid, 'pr_ens_opt_f', 'double', [dimid2 dimid0]);
varid13  = netcdf.defVar(myncid, 'pr_ens_opt_xi', 'double', [dimid2 dimid0]);
varid21  = netcdf.defVar(myncid, 'tasmax_ens_opt_stat', 'double', [dimid1 dimid0]);
varid22  = netcdf.defVar(myncid, 'tasmax_ens_opt_f', 'double', [dimid2 dimid0]);
varid23  = netcdf.defVar(myncid, 'tasmax_ens_opt_xi', 'double', [dimid2 dimid0]);
varid31  = netcdf.defVar(myncid, 'spna_obs', 'double', [dimid4]);
varid32  = netcdf.defVar(myncid, 'spna_ens_opt', 'double', [dimid5 dimid3]);
varid33  = netcdf.defVar(myncid, 'spna_ens_opt_mme', 'double', [dimid5]);
varid4   = netcdf.defVar(myncid, 'ba_frac_prob', 'double', [dimid6]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid11, pr_ens_opt_stat);
netcdf.putVar(myncid, varid12, pr_ens_opt_f);
netcdf.putVar(myncid, varid13, pr_ens_opt_xi);
netcdf.putVar(myncid, varid21, tasmax_ens_opt_stat);
netcdf.putVar(myncid, varid22, tasmax_ens_opt_f);
netcdf.putVar(myncid, varid23, tasmax_ens_opt_xi);
netcdf.putVar(myncid, varid31, spna(42:62));
netcdf.putVar(myncid, varid32, spna_ens_opt(1:50,:));
netcdf.putVar(myncid, varid33, mean(spna_ens_opt(1:50,:),2));
netcdf.putVar(myncid, varid4,  ba_frac_prob);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lon = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','lon');
lat = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','lat');
pre = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.pre.dat.nc','pre');
tmx = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\cru_ts4.07.1901.2022.tmx.dat.nc','tmx');

pre = pre(:,:,1201:1452);
tmx = tmx(:,:,1201:1452);

pre = reshape(pre,[720*360 252]);
tmx = reshape(tmx,[720*360 252]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);

pre = pre(lct,:);
tmx = tmx(lct,:);
lon1d_clim = lon1d(lct);
lat1d_clim = lat1d(lct);

pre  = reshape(pre,[1500 12 21]);
tmx  = reshape(tmx,[1500 12 21]);

pre_ann = squeeze(sum(pre(:,4:9,:),2));
tmx_ann = squeeze(mean(tmx(:,4:9,:),2));

clear lat lat1d lat2d lon lon1d lon2d lct pre tmx 

lon = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lon');
lat = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','lat');
yr  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','year');
ba  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','burn_area');
ta  = ncread('H:\WORKS\29-Siberian_Fire\ext_fig1\MCD64A1_burned_area_ne.nc','total_area');

ba  = 100*ba./ta;
ba  = ba(:,:,:,1:21);
ba  = reshape(ba,[144*16 12*21]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[144*16 1]);
lat1d  = reshape(lat2d,[144*16 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);

ba = ba(lct,:);
lon1d_ba = lon1d(lct);
lat1d_ba = lat1d(lct);

ba = reshape(ba,[60 12 21]);
ba_ann = squeeze(sum(ba(:,4:9,:),2));

clear lon lon1d lon2d lat lat1d lat2d yr ba ta

ba_ann = reshape(ba_ann,[12 5 21]);
lat = reshape(lat1d_ba,[12 5]);
lon = reshape(lon1d_ba,[12 5]);

for j = 1:size(lon,2)
    for i = 1:size(lon,1)
        lct = find(lon1d_clim>=lon(i,j)-1.25&lon1d_clim<lon(i,j)+1.25&...
                   lat1d_clim>=lat(i,j)-1.25&lat1d_clim<lat(i,j)+1.25);
        pre_ann_25deg(i,j,:) = mean(pre_ann(lct,:),1);
        tmx_ann_25deg(i,j,:) = mean(tmx_ann(lct,:),1);
        
        clear lct
    end
end

for j = 1:size(lon,2)
    for i = 1:size(lon,1)
    pre_ann_sigma(i,j,:) = (pre_ann_25deg(i,j,:)-mean(pre_ann_25deg(i,j,:)))/std(pre_ann_25deg(i,j,:));
    tmx_ann_sigma(i,j,:) = (tmx_ann_25deg(i,j,:)-mean(tmx_ann_25deg(i,j,:)))/std(tmx_ann_25deg(i,j,:));
    end
end        

pre_ann_diff = diff(pre_ann_sigma,1,3);
tmx_ann_diff = diff(tmx_ann_sigma,1,3);
ba_ann_diff  = diff(ba_ann,1,3);

for j = 1:size(lon,2)
    for i = 1:size(lon,1)
    b = polyfit(squeeze(pre_ann_diff(i,j,:)),squeeze(ba_ann_diff(i,j,:)),1);
    pre_sens(i,j) = b(1);
    clear b

    b = polyfit(squeeze(tmx_ann_diff(i,j,:)),squeeze(ba_ann_diff(i,j,:)),1);
    tmx_sens(i,j) = b(1);
    clear b
    end
end

myncid  = netcdf.create('figs1.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 12);
dimid2  = netcdf.defDim(myncid, 'lat', 5);
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pre_sens', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'tmx_sens', 'double', [dimid1 dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon(:,1));
netcdf.putVar(myncid, varid2, lat(1,:)');
netcdf.putVar(myncid, varid3, pre_sens);
netcdf.putVar(myncid, varid4, tmx_sens);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

pr_sib_trd(1)      = mean(ncread('fig3a-accessesm.nc','pr_sib_trd'));
pr_sib_trd(2)      = mean(ncread('fig3a-canesm2.nc','pr_sib_trd'));
pr_sib_trd(3)      = mean(ncread('fig3a-canesm5.nc','pr_sib_trd'));
pr_sib_trd(4)      = mean(ncread('fig3a-miroc6.nc','pr_sib_trd'));
pr_sib_trd(5)      = mean(ncread('fig3a-mpi-ge.nc','pr_sib_trd'));
pr_sib_trd(6)      = mean(ncread('fig3a-mpiesm12.nc','pr_sib_trd'));

pr_sib_trd_p5(1)   = prctile(ncread('fig3a-accessesm.nc','pr_sib_trd'),5);
pr_sib_trd_p5(2)   = prctile(ncread('fig3a-canesm2.nc','pr_sib_trd'),5);
pr_sib_trd_p5(3)   = prctile(ncread('fig3a-canesm5.nc','pr_sib_trd'),5);
pr_sib_trd_p5(4)   = prctile(ncread('fig3a-miroc6.nc','pr_sib_trd'),5);
pr_sib_trd_p5(5)   = prctile(ncread('fig3a-mpi-ge.nc','pr_sib_trd'),5);
pr_sib_trd_p5(6)   = prctile(ncread('fig3a-mpiesm12.nc','pr_sib_trd'),5);

pr_sib_trd_p95(1)  = prctile(ncread('fig3a-accessesm.nc','pr_sib_trd'),95);
pr_sib_trd_p95(2)  = prctile(ncread('fig3a-canesm2.nc','pr_sib_trd'),95);
pr_sib_trd_p95(3)  = prctile(ncread('fig3a-canesm5.nc','pr_sib_trd'),95);
pr_sib_trd_p95(4)  = prctile(ncread('fig3a-miroc6.nc','pr_sib_trd'),95);
pr_sib_trd_p95(5)  = prctile(ncread('fig3a-mpi-ge.nc','pr_sib_trd'),95);
pr_sib_trd_p95(6)  = prctile(ncread('fig3a-mpiesm12.nc','pr_sib_trd'),95);

pr_sib_trd_sigma(1)      = mean(ncread('fig3a-accessesm.nc','pr_sib_trd_sigma'));
pr_sib_trd_sigma(2)      = mean(ncread('fig3a-canesm2.nc','pr_sib_trd_sigma'));
pr_sib_trd_sigma(3)      = mean(ncread('fig3a-canesm5.nc','pr_sib_trd_sigma'));
pr_sib_trd_sigma(4)      = mean(ncread('fig3a-miroc6.nc','pr_sib_trd_sigma'));
pr_sib_trd_sigma(5)      = mean(ncread('fig3a-mpi-ge.nc','pr_sib_trd_sigma'));
pr_sib_trd_sigma(6)      = mean(ncread('fig3a-mpiesm12.nc','pr_sib_trd_sigma'));

pr_sib_trd_sigma_p5(1)   = prctile(ncread('fig3a-accessesm.nc','pr_sib_trd_sigma'),5);
pr_sib_trd_sigma_p5(2)   = prctile(ncread('fig3a-canesm2.nc','pr_sib_trd_sigma'),5);
pr_sib_trd_sigma_p5(3)   = prctile(ncread('fig3a-canesm5.nc','pr_sib_trd_sigma'),5);
pr_sib_trd_sigma_p5(4)   = prctile(ncread('fig3a-miroc6.nc','pr_sib_trd_sigma'),5);
pr_sib_trd_sigma_p5(5)   = prctile(ncread('fig3a-mpi-ge.nc','pr_sib_trd_sigma'),5);
pr_sib_trd_sigma_p5(6)   = prctile(ncread('fig3a-mpiesm12.nc','pr_sib_trd_sigma'),5);

pr_sib_trd_sigma_p95(1)  = prctile(ncread('fig3a-accessesm.nc','pr_sib_trd_sigma'),95);
pr_sib_trd_sigma_p95(2)  = prctile(ncread('fig3a-canesm2.nc','pr_sib_trd_sigma'),95);
pr_sib_trd_sigma_p95(3)  = prctile(ncread('fig3a-canesm5.nc','pr_sib_trd_sigma'),95);
pr_sib_trd_sigma_p95(4)  = prctile(ncread('fig3a-miroc6.nc','pr_sib_trd_sigma'),95);
pr_sib_trd_sigma_p95(5)  = prctile(ncread('fig3a-mpi-ge.nc','pr_sib_trd_sigma'),95);
pr_sib_trd_sigma_p95(6)  = prctile(ncread('fig3a-mpiesm12.nc','pr_sib_trd_sigma'),95);

pr_nor_trd(1)      = mean(ncread('fig3a-accessesm.nc','pr_nor_trd'));
pr_nor_trd(2)      = mean(ncread('fig3a-canesm2.nc','pr_nor_trd'));
pr_nor_trd(3)      = mean(ncread('fig3a-canesm5.nc','pr_nor_trd'));
pr_nor_trd(4)      = mean(ncread('fig3a-miroc6.nc','pr_nor_trd'));
pr_nor_trd(5)      = mean(ncread('fig3a-mpi-ge.nc','pr_nor_trd'));
pr_nor_trd(6)      = mean(ncread('fig3a-mpiesm12.nc','pr_nor_trd'));

pr_nor_trd_p5(1)   = prctile(ncread('fig3a-accessesm.nc','pr_nor_trd'),5);
pr_nor_trd_p5(2)   = prctile(ncread('fig3a-canesm2.nc','pr_nor_trd'),5);
pr_nor_trd_p5(3)   = prctile(ncread('fig3a-canesm5.nc','pr_nor_trd'),5);
pr_nor_trd_p5(4)   = prctile(ncread('fig3a-miroc6.nc','pr_nor_trd'),5);
pr_nor_trd_p5(5)   = prctile(ncread('fig3a-mpi-ge.nc','pr_nor_trd'),5);
pr_nor_trd_p5(6)   = prctile(ncread('fig3a-mpiesm12.nc','pr_nor_trd'),5);

pr_nor_trd_p95(1)  = prctile(ncread('fig3a-accessesm.nc','pr_nor_trd'),95);
pr_nor_trd_p95(2)  = prctile(ncread('fig3a-canesm2.nc','pr_nor_trd'),95);
pr_nor_trd_p95(3)  = prctile(ncread('fig3a-canesm5.nc','pr_nor_trd'),95);
pr_nor_trd_p95(4)  = prctile(ncread('fig3a-miroc6.nc','pr_nor_trd'),95);
pr_nor_trd_p95(5)  = prctile(ncread('fig3a-mpi-ge.nc','pr_nor_trd'),95);
pr_nor_trd_p95(6)  = prctile(ncread('fig3a-mpiesm12.nc','pr_nor_trd'),95);

pr_nor_trd_sigma(1)      = mean(ncread('fig3a-accessesm.nc','pr_nor_trd_sigma'));
pr_nor_trd_sigma(2)      = mean(ncread('fig3a-canesm2.nc','pr_nor_trd_sigma'));
pr_nor_trd_sigma(3)      = mean(ncread('fig3a-canesm5.nc','pr_nor_trd_sigma'));
pr_nor_trd_sigma(4)      = mean(ncread('fig3a-miroc6.nc','pr_nor_trd_sigma'));
pr_nor_trd_sigma(5)      = mean(ncread('fig3a-mpi-ge.nc','pr_nor_trd_sigma'));
pr_nor_trd_sigma(6)      = mean(ncread('fig3a-mpiesm12.nc','pr_nor_trd_sigma'));

pr_nor_trd_sigma_p5(1)   = prctile(ncread('fig3a-accessesm.nc','pr_nor_trd_sigma'),5);
pr_nor_trd_sigma_p5(2)   = prctile(ncread('fig3a-canesm2.nc','pr_nor_trd_sigma'),5);
pr_nor_trd_sigma_p5(3)   = prctile(ncread('fig3a-canesm5.nc','pr_nor_trd_sigma'),5);
pr_nor_trd_sigma_p5(4)   = prctile(ncread('fig3a-miroc6.nc','pr_nor_trd_sigma'),5);
pr_nor_trd_sigma_p5(5)   = prctile(ncread('fig3a-mpi-ge.nc','pr_nor_trd_sigma'),5);
pr_nor_trd_sigma_p5(6)   = prctile(ncread('fig3a-mpiesm12.nc','pr_nor_trd_sigma'),5);

pr_nor_trd_sigma_p95(1)  = prctile(ncread('fig3a-accessesm.nc','pr_nor_trd_sigma'),95);
pr_nor_trd_sigma_p95(2)  = prctile(ncread('fig3a-canesm2.nc','pr_nor_trd_sigma'),95);
pr_nor_trd_sigma_p95(3)  = prctile(ncread('fig3a-canesm5.nc','pr_nor_trd_sigma'),95);
pr_nor_trd_sigma_p95(4)  = prctile(ncread('fig3a-miroc6.nc','pr_nor_trd_sigma'),95);
pr_nor_trd_sigma_p95(5)  = prctile(ncread('fig3a-mpi-ge.nc','pr_nor_trd_sigma'),95);
pr_nor_trd_sigma_p95(6)  = prctile(ncread('fig3a-mpiesm12.nc','pr_nor_trd_sigma'),95);

myncid  = netcdf.create('figs2.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'model', 6);
varid11  = netcdf.defVar(myncid, 'pr_sib_trd', 'double', dimid1);
varid12  = netcdf.defVar(myncid, 'pr_nor_trd', 'double', dimid1);
varid13  = netcdf.defVar(myncid, 'pr_sib_trd_p5', 'double', dimid1);
varid14  = netcdf.defVar(myncid, 'pr_nor_trd_p5', 'double', dimid1);
varid15  = netcdf.defVar(myncid, 'pr_sib_trd_p95', 'double', dimid1);
varid16  = netcdf.defVar(myncid, 'pr_nor_trd_p95', 'double', dimid1);
varid21  = netcdf.defVar(myncid, 'pr_sib_trd_sigma', 'double', dimid1);
varid22  = netcdf.defVar(myncid, 'pr_nor_trd_sigma', 'double', dimid1);
varid23  = netcdf.defVar(myncid, 'pr_sib_trd_sigma_p5', 'double', dimid1);
varid24  = netcdf.defVar(myncid, 'pr_nor_trd_sigma_p5', 'double', dimid1);
varid25  = netcdf.defVar(myncid, 'pr_sib_trd_sigma_p95', 'double', dimid1);
varid26  = netcdf.defVar(myncid, 'pr_nor_trd_sigma_p95', 'double', dimid1);
netcdf.endDef(myncid);

netcdf.putVar(myncid, varid11, pr_sib_trd);
netcdf.putVar(myncid, varid12, pr_nor_trd);
netcdf.putVar(myncid, varid13, pr_sib_trd_p5);
netcdf.putVar(myncid, varid14, pr_nor_trd_p5);
netcdf.putVar(myncid, varid15, pr_sib_trd_p95);
netcdf.putVar(myncid, varid16, pr_nor_trd_p95);
netcdf.putVar(myncid, varid21, pr_sib_trd_sigma);
netcdf.putVar(myncid, varid22, pr_nor_trd_sigma);
netcdf.putVar(myncid, varid23, pr_sib_trd_sigma_p5);
netcdf.putVar(myncid, varid24, pr_nor_trd_sigma_p5);
netcdf.putVar(myncid, varid25, pr_sib_trd_sigma_p95);
netcdf.putVar(myncid, varid26, pr_nor_trd_sigma_p95);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lat    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','lat');
lon    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','lon');
time   = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','time');
pre    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','precip'); % 1891-2019
pre = pre(:,:,1321:1548); % 2001-2019

time2  = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.first.mon.total.1x1.nc','time');
pre2   = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.first.mon.total.1x1.nc','precip'); % 2012-2023
pre2 = pre2(:,:,97:120); % 2020-2021

pre = cat(3,pre,pre2);
pre(pre<0) = NaN;

pre = reshape(pre,[360*180 252]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[360*180 1]);
lat1d  = reshape(lat2d,[360*180 1]);

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');

for j=1:length(lat)
for i=1:length(lon)
    lct = find(lat1d>=lat(j)-1&lat1d<lat(j)+1&lon1d>=lon(i)-1&lon1d<lon(i)+1);
    pre_2deg(i,j,:) = mean(pre(lct,:),1);
end
disp(j)
end

pre = reshape(pre_2deg,[180 90 12*21]); % 1960-2021

nansum = sum(isnan(pre),3);

pre = reshape(pre,[180 90 12 21]);
pre_ann = squeeze(sum(pre(:,:,4:9,:),3));

year  = [2001:2021]';

for j = 1:size(pre_ann,2)
for i = 1:size(pre_ann,1)
      if nansum(i,j) ==0
      [b,bint,r,rint,stats] = regress(squeeze(pre_ann(i,j,:)),[ones(size(year)) year]);
      pre_trd(i,j) = b(2)*10;
      pre_trd_pval(i,j) = stats(3);
      else
      pre_trd(i,j) = -999;
      pre_trd_pval(i,j) = 1;
      end
end
end

pre_trd_pval = reshape(pre_trd_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

lonsig = lon1d(find(pre_trd_pval<=0.05));
latsig = lat1d(find(pre_trd_pval<=0.05));

myncid  = netcdf.create('figs3a.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));

varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'pre_trd', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);

netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, pre_trd);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

diri = 'H:\WORKS\29-Siberian_Fire\figs3\';
plev   = ncread(strcat(diri,'hgt.mon.mean.nc'),'level');
zg   = ncread(strcat(diri,'hgt.mon.mean.nc'),'hgt'); %1979-2023
zg   = squeeze(zg(:,:,6,1+12*22:12*43)); % 2001-2021

lat  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lat');
lon  = ncread(strcat(diri,'hgt.mon.mean.nc'),'lon');
[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for k=1:size(zg,3)
    zg_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),zg(:,end:-1:1,k),Xi,Yi);
    disp(k)
end
clear zg

zg_2deg(end,:,:) = zg_2deg(end-1,:,:);

zg = reshape(zg_2deg,[180 90 12 21]); % 2001-2021
zg = zg - repmat(mean(zg,1),[180 1 1 1]); % eddy geopotential height

zg_ann = squeeze(mean(zg(:,:,4:9,:),3));% 2001-2021

year  = [2001:2021]';

for j = 1:size(zg_ann,2)
for i = 1:size(zg_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(zg_ann(i,j,:)),[ones(size(year)) year]);
      zg_trd(i,j) = b(2)*10;
      zg_trd_pval(i,j) = stats(3);
end
end

zg_trd_pval = reshape(zg_trd_pval,[180*90 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*90 1]);
lat1d  = reshape(lat2d,[180*90 1]);

lonsig = lon1d(find(zg_trd_pval<=0.05));
latsig = lat1d(find(zg_trd_pval<=0.05));

myncid  = netcdf.create('figs3b.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
dimid3  = netcdf.defDim(myncid, 'sig', length(lonsig));

varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'zg_trd', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'lonsig', 'double', [dimid3]);
varid5  = netcdf.defVar(myncid, 'latsig', 'double', [dimid3]);
netcdf.endDef(myncid);

netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, zg_trd);
netcdf.putVar(myncid, varid4, lonsig);
netcdf.putVar(myncid, varid5, latsig);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lat    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','lat');
lon    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','lon');
time   = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','time');
pre    = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.mon.total.1x1.v2020.nc','precip'); % 1891-2019
pre  = pre(:,:,829:1548); % 1960-2019

pre2   = ncread('H:\WORKS\29-Siberian_Fire\figs3\precip.first.mon.total.1x1.nc','precip'); % 2012-2023
pre2 = pre2(:,:,97:120); % 2020-2021

pre  = cat(3,pre,pre2);
pre(pre<0) = NaN;

pre = reshape(pre,[360*180 744]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[360 1]);
weight = reshape(weight,[360*180 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[360*180 1]);
lat1d  = reshape(lat2d,[360*180 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);
pre_avg = mean(pre.*repmat(weights,[1 744]),1);
pre_avg = reshape(pre_avg,[12 62]);
pre_ann = sum(pre_avg(4:9,:),1);
pre_ann = pre_ann - mean(pre_ann);

sst  = ncread('H:\WORKS\29-Siberian_Fire\figs4\HadSST.4.0.1.0_median.nc','tos'); % since 1850
lat  = ncread('H:\WORKS\29-Siberian_Fire\figs4\HadSST.4.0.1.0_median.nc','latitude'); 
lon  = ncread('H:\WORKS\29-Siberian_Fire\figs4\HadSST.4.0.1.0_median.nc','longitude');

sst  = sst(:,:,1321:2064); % 1960-2021
sst  = reshape(sst,[size(sst,1)*size(sst,2) size(sst,3)]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[72*36 1]);
lat1d  = reshape(lat2d,[72*36 1]);

weight = cos(pi*lat/180.0);
weight = weight/mean(weight);
weight = repmat(weight',[length(lon) 1]);
weight = reshape(weight,[size(weight,1)*size(weight,2) 1]);

for i = 1:size(sst,2)
   lct        = find(isnan(sst(:,i))==0&lat1d>=-45&lat1d<=65);
   weight2    = weight(lct)/mean(weight(lct));
   sst_glb(i) = mean(sst(lct,i).*weight2);
   clear lct weight2
end
clear sst weight

sst  = ncread('H:\WORKS\29-Siberian_Fire\figs4\HadSST.4.0.1.0_median.nc','tos'); % since 1850
sst  = sst(27:34,28:31,1321:2064); % 50-10W; 45-65N; 1960-2021
sst  = reshape(sst,[size(sst,1)*size(sst,2) size(sst,3)]);
nan_sum = sum(isnan(sst),2);
lct     = find(nan_sum~=size(sst,2));

sst  = sst(lct,:);
sst_res = zeros(size(sst));

for i = 1:size(sst,1)
   lct2 = find(isnan(sst(i,:))==0);
   b = regress(sst(i,lct2)',[ones(size(lct2))' sst_glb(lct2)']);
   sst_res(i,:) = sst(i,:) - (b(2)*sst_glb+b(1));
   clear lct2 b
end

weight = cos(pi*lat(28:31)./180.0);
weight = weight/mean(weight);
weight = repmat(weight',[8 1]);
weight = reshape(weight,[size(weight,1)*size(weight,2) 1]);
weight2 = weight(lct)/mean(weight(lct));

for i = 1:size(sst_res,2)
   lct3 = find(isnan(sst_res(:,i))==0);
   weight3 = weight2(lct3)/mean(weight2(lct3));
   spna(i) = mean(sst_res(lct3,i).*weight3);
   clear lct3 weight3
end

spna = reshape(spna,[12 62]);
spna = mean(spna(4:9,:),1);
spna = movmean(spna,9);

[R,P] = corrcoef(spna,pre_ann)

myncid  = netcdf.create('figs4.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'yr', 62);
varid1  = netcdf.defVar(myncid, 'spna', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_ann', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pre_ann_movavg', 'double', [dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, spna);
netcdf.putVar(myncid, varid2, pre_ann);
netcdf.putVar(myncid, varid3, movmean(pre_ann,9));
netcdf.close(myncid);
%**************************************************************************
clear all
clc

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

diri = 'H:\WORKS\29-Siberian_Fire\figure4\';
lat   = ncread(strcat(diri,'uwnd.mon.mean.nc'),'lat');
lon   = ncread(strcat(diri,'uwnd.mon.mean.nc'),'lon');
plev  = ncread(strcat(diri,'uwnd.mon.mean.nc'),'level');

ua   = ncread(strcat(diri,'uwnd.mon.mean.nc'),'uwnd');
ua   = ua(:,:,1:10,145:12*74); % 1960-2021
wa   = ncread(strcat(diri,'omega.mon.mean.nc'),'omega');
wa   = wa(:,:,1:10,145:12*74); % 1960-2021

[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi]  = meshgrid(lat',lon');

for k=1:size(ua,4)
    for p=1:size(ua,3)
    ua_2deg(:,:,p,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),ua(:,end:-1:1,p,k),Xi,Yi);
    wa_2deg(:,:,p,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),wa(:,end:-1:1,p,k),Xi,Yi);
    end
    disp(k)
end

ua_2deg(end,:,:) = ua_2deg(end-1,:,:);
ua = reshape(ua_2deg,[180 90 10 12 62]); % 1960-2021
ua_ann = squeeze(mean(ua(:,:,:,4:9,:),4));% 1960-2021

wa_2deg(end,:,:) = wa_2deg(end-1,:,:);
wa = reshape(wa_2deg,[180 90 10 12 62]); % 1960-2021
wa_ann = squeeze(mean(wa(:,:,:,4:9,:),4));% 1960-2021

for p = 1:size(ua_ann,3)
for j = 1:size(ua_ann,2)
for i = 1:size(ua_ann,1)
      [b,bint,r,rint,stats] = regress(squeeze(ua_ann(i,j,p,:)),[ones(size(spna)) spna]);
      ua_reg_spna(i,j,p) = b(2);

      [b,bint,r,rint,stats] = regress(squeeze(wa_ann(i,j,p,:)),[ones(size(spna)) spna]);
      wa_reg_spna(i,j,p) = b(2);
end
end
disp(p)
end

ua_reg_spna = squeeze(mean(cat(1,ua_reg_spna(166:180,73:79,:),ua_reg_spna(1:90,73:79,:)),2));
wa_reg_spna = squeeze(mean(cat(1,wa_reg_spna(166:180,73:79,:),wa_reg_spna(1:90,73:79,:)),2));

contourf(wa_reg_spna)

myncid  = netcdf.create('figs5.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 105);
dimid2  = netcdf.defDim(myncid, 'plev', 10);

varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'plev', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'ua_reg_spna', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'wa_reg_spna', 'double', [dimid1 dimid2]);
netcdf.endDef(myncid);

netcdf.putVar(myncid, varid1, [lon(166:180)'-360 lon(1:90)']');
netcdf.putVar(myncid, varid2, plev(1:10));
netcdf.putVar(myncid, varid3, ua_reg_spna);
netcdf.putVar(myncid, varid4, wa_reg_spna);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lon  = ncread('H:\WORKS\29-Siberian_Fire\figs6\TN2001-Fx.monthly.nc','lon');
lat  = ncread('H:\WORKS\29-Siberian_Fire\figs6\TN2001-Fx.monthly.nc','lat');
plev = ncread('H:\WORKS\29-Siberian_Fire\figs6\TN2001-Fx.monthly.nc','level');
Fx   = ncread('H:\WORKS\29-Siberian_Fire\figs6\TN2001-Fx.monthly.nc','Fx');
Fy   = ncread('H:\WORKS\29-Siberian_Fire\figs6\TN2001-Fy.monthly.nc','Fy');
psi  = ncread('H:\WORKS\29-Siberian_Fire\figs6\psidev.monthly.nc','psidev');

Fx   = squeeze(Fx(:,:,8,:));  % 1960-2021
Fy   = squeeze(Fy(:,:,8,:));  % 1960-2021
psi  = squeeze(psi(:,:,8,:)); % 1960-2021

[X,Y]  = meshgrid(lat',lon');
clear lat lon

lat     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lat');
lon     = ncread('H:\WORKS\29-Siberian_Fire\figure2\land_sea_mask.nc','lon');
[Xi,Yi] = meshgrid(lat',lon');

for k=1:size(Fx,3)
    Fx_2deg(:,:,k)  = interp2(X(:,end:-1:1),Y(:,end:-1:1),Fx(:,end:-1:1,k),Xi,Yi);
    Fy_2deg(:,:,k)  = interp2(X(:,end:-1:1),Y(:,end:-1:1),Fy(:,end:-1:1,k),Xi,Yi);
    psi_2deg(:,:,k) = interp2(X(:,end:-1:1),Y(:,end:-1:1),psi(:,end:-1:1,k),Xi,Yi);
    disp(k)
end

clear Fx Fy psi

Fx_2deg(end,:,:)  = Fx_2deg(end-1,:,:);
Fy_2deg(end,:,:)  = Fy_2deg(end-1,:,:);
psi_2deg(end,:,:) = psi_2deg(end-1,:,:);

Fx_ann  = squeeze(mean(Fx_2deg(:,:,4:9),3)); % 1960-2021
Fy_ann  = squeeze(mean(Fy_2deg(:,:,4:9),3)); % 1960-2021
psi_ann = squeeze(mean(psi_2deg(:,:,4:9),3));% 1960-2021

myncid  = netcdf.create('figs6.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 180);
dimid2  = netcdf.defDim(myncid, 'lat', 90);
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'Fx_ann', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'Fy_ann', 'double', [dimid1 dimid2]);
varid5  = netcdf.defVar(myncid, 'psi_ann', 'double', [dimid1 dimid2]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, Fx_ann);
netcdf.putVar(myncid, varid4, Fy_ann);
netcdf.putVar(myncid, varid5, psi_ann);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

lat = ncread('H:\WORKS\29-Siberian_Fire\figs7\response.spna.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figs7\response.spna.nc','lon');
lev = ncread('H:\WORKS\29-Siberian_Fire\figs7\response.spna.nc','lev');
[X,Y] = meshgrid(lat,lon);

u   = ncread('H:\WORKS\29-Siberian_Fire\figs7\response.spna.nc','u');
v   = ncread('H:\WORKS\29-Siberian_Fire\figs7\response.spna.nc','v');

u700 = mean(u(:,:,5,21:50),4);
v700 = mean(v(:,:,5,21:50),4);

u700(find(isnan(u700)==1)) = -999;
v700(find(isnan(v700)==1)) = -999;

lat0 = 55;
lon0 = 330;
semimajor = 20;
ecc = axes2ecc(20,10);
[late,lone] = ellipse1(lat0,lon0,[semimajor ecc],90);

myncid  = netcdf.create('lbm_response.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'lon', 128);
dimid2  = netcdf.defDim(myncid, 'lat', 64);
dimid3  = netcdf.defDim(myncid, 'ell', 100);
varid1  = netcdf.defVar(myncid, 'lon', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'lat', 'double', [dimid2]);
varid3  = netcdf.defVar(myncid, 'u700', 'double', [dimid1 dimid2]);
varid4  = netcdf.defVar(myncid, 'v700', 'double', [dimid1 dimid2]);
varid5  = netcdf.defVar(myncid, 'lone', 'double', [dimid3]);
varid6  = netcdf.defVar(myncid, 'late', 'double', [dimid3]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, lon);
netcdf.putVar(myncid, varid2, lat);
netcdf.putVar(myncid, varid3, u700);
netcdf.putVar(myncid, varid4, v700);
netcdf.putVar(myncid, varid5, lone);
netcdf.putVar(myncid, varid6, late);
netcdf.close(myncid);
%**************************************************************************
clear all
clc

load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\spna_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\pr_sib_fgoals_g3.mat
load H:\WORKS\29-Siberian_Fire\figure2\FGOALS-g3\tasmax_sib_fgoals_g3.mat

[lcti lctj] = find(tasmax_ann_sib>1000);
tasmax_ann_sib(lcti,lctj) = (tasmax_ann_sib(lcti-1,lctj)+tasmax_ann_sib(lcti+1,lctj))*0.5; 

spna_ens = spna(101:206,:); % 1950-2055

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 110]);
pr_sib_ens     = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]); % 1950-2055
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = pr_sib_forced;

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 110]);
tasmax_sib_ens     = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = tasmax_sib_forced;

clear spna pr_ann_sib pr_ann_sib_int pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\spna_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\pr_sib_accessesm.mat
load H:\WORKS\29-Siberian_Fire\figure2\ACCESS-ESM1-5\tasmax_sib_accessesm.mat

spna_ens = cat(2,spna_ens,spna(101:206,:)); % 1950-2055

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 40]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 40]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\spna_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\pr_sib_canesm5.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM5\tasmax_sib_canesm5.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\spna_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\pr_sib_miroc6.mat
load H:\WORKS\29-Siberian_Fire\figure2\MIROC6\tasmax_sib_miroc6.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\spna_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\pr_sib_mpiesm12.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-ESM1-2-LR\tasmax_sib_mpiesm12.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 44]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 44]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\spna_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\pr_sib_canesm2.mat
load H:\WORKS\29-Siberian_Fire\figure2\CanESM2\tasmax_sib_canesm2.mat

spna_ens = cat(2,spna_ens,spna(1:106,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 50]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tasmax_ann_sib - repmat(mean(tasmax_ann_sib,2),[1 50]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tasmax_ann_sib(1:106,:)-repmat(mean(tasmax_ann_sib(52:72,:),1),[106 1]))./repmat(std(tasmax_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\spna_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\pr_sib_mpi_ge.mat
load H:\WORKS\29-Siberian_Fire\figure2\MPI-GE\tas_sib_mpi_ge.mat

spna_ens = cat(2,spna_ens,spna(101:206,:));

pr_ann_sib_int = pr_ann_sib - repmat(mean(pr_ann_sib,2),[1 100]);
pr_ann_sib_sig = (pr_ann_sib_int(1:106,:)-repmat(mean(pr_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]);
pr_sib_ens = cat(2,pr_sib_ens,pr_ann_sib_sig);
pr_sib_forced  = mean((pr_ann_sib(1:106,:)-repmat(mean(pr_ann_sib(52:72,:),1),[106 1]))./repmat(std(pr_ann_sib(52:72,:),0,1),[106 1]),2);
pr_sib_ens2    = cat(2,pr_sib_ens2,pr_sib_forced);

tasmax_ann_sib_int = tas_ann_sib - repmat(mean(tas_ann_sib,2),[1 100]);
tasmax_ann_sib_sig = (tasmax_ann_sib_int(1:106,:)-repmat(mean(tasmax_ann_sib_int(52:72,:),1),[106 1]))./repmat(std(tas_ann_sib(52:72,:),0,1),[106 1]);
tasmax_sib_ens = cat(2,tasmax_sib_ens,tasmax_ann_sib_sig);
tasmax_sib_forced  = mean((tas_ann_sib(1:106,:)-repmat(mean(tas_ann_sib(52:72,:),1),[106 1]))./repmat(std(tas_ann_sib(52:72,:),0,1),[106 1]),2);
tasmax_sib_ens2    = cat(2,tasmax_sib_ens2,tasmax_sib_forced);

clear spna pr_ann_sib pr_ann_sib_int pr_ann_sib_sig pr_sib_forced tasmax_ann_sib tasmax_ann_sib_int tasmax_ann_sib_sig tasmax_sib_forced ts_ann_res_spna ts_ann_spna

spna = ncread('H:\WORKS\29-Siberian_Fire\figure3\fig3e.nc','spna');

lct = ncread('H:\WORKS\29-Siberian_Fire\figure5\fig5a.nc','lct');
spna_ens = spna_ens(:,lct);

b0 = polyfit([2001:2021]',spna(42:62),1);
spna_trd = b0(1)*10;

for i = 1:length(lct)
for j = 1:52
    b1 = polyfit([2021:2050]',spna_ens(j+20:j+49,i),1);
    spna_ens_trd(j) = b1(1)*10;
    clear b1
end

lct_max_corr         = find(spna_ens_trd==max(spna_ens_trd));
spna_ens_opt(:,i)     = spna_ens(lct_max_corr:lct_max_corr+54,i);
pr_ens_opt(:,i)      = pr_sib_ens(lct_max_corr:lct_max_corr+54,i);
tasmax_ens_opt(:,i)  = tasmax_sib_ens(lct_max_corr:lct_max_corr+54,i);

clear lct_max_corr
end

pr_sib_ens2 = mean(pr_sib_ens2,2);
tasmax_sib_ens2 = mean(tasmax_sib_ens2,2);

pr_sib_ens2_hist     = pr_sib_ens2(11:72);
tasmax_sib_ens2_hist = tasmax_sib_ens2(11:72);

pre_ann    = ncread('H:\WORKS\29-Siberian_Fire\figure5\obs_pr_tasmax.nc','pre_ann');
tasmax_ann = ncread('H:\WORKS\29-Siberian_Fire\figure5\obs_pr_tasmax.nc','tmx_ann');

pre_ann    = (pre_ann-mean(pre_ann(42:62)))/std(pre_ann(42:62),0,1);
tasmax_ann = (tasmax_ann-mean(tasmax_ann(42:62)))/std(tasmax_ann(42:62),0,1);

scale_pr     = polyfit(pr_sib_ens2_hist,pre_ann,1);
scale_tasmax = polyfit(tasmax_sib_ens2_hist,tasmax_ann,1);

pr_sib_ens = pr_sib_ens2(52:106)*scale_pr(1); % 2001-2050
tasmax_sib_ens = tasmax_sib_ens2(52:106)*scale_tasmax(1);


polyfit([2021:2040],mean(pr_ens_opt(21:40,:),2),1)
polyfit([2021:2050],pr_sib_ens(21:50),1)

polyfit([2021:2050],mean(tasmax_ens_opt(21:50,:),2),1)
polyfit([2021:2050],tasmax_sib_ens(21:50),1)

myncid  = netcdf.create('figs8.nc', 'NC_NOCLOBBER');
dimid0  = netcdf.defDim(myncid, 'year', 30);
varid1  = netcdf.defVar(myncid, 'pr_sib_ens', 'double', [dimid0]);
varid2  = netcdf.defVar(myncid, 'tasmax_sib_ens', 'double', [dimid0]);
varid11  = netcdf.defVar(myncid, 'pr_ens_opt', 'double', [dimid0]);
varid12  = netcdf.defVar(myncid, 'pr_ens_opt_p5', 'double', [dimid0]);
varid13  = netcdf.defVar(myncid, 'pr_ens_opt_p95', 'double', [dimid0]);
varid21  = netcdf.defVar(myncid, 'tasmax_ens_opt', 'double', [dimid0]);
varid22  = netcdf.defVar(myncid, 'tasmax_ens_opt_p5', 'double', [dimid0]);
varid23  = netcdf.defVar(myncid, 'tasmax_ens_opt_p95', 'double', [dimid0]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, pr_sib_ens(21:50));
netcdf.putVar(myncid, varid2, tasmax_sib_ens(21:50));
netcdf.putVar(myncid, varid11, mean(pr_ens_opt(21:50,:),2));
netcdf.putVar(myncid, varid12, prctile(pr_ens_opt(21:50,:),5,2));
netcdf.putVar(myncid, varid13, prctile(pr_ens_opt(21:50,:),95,2));
netcdf.putVar(myncid, varid21, mean(tasmax_ens_opt(21:50,:),2));
netcdf.putVar(myncid, varid22, prctile(tasmax_ens_opt(21:50,:),5,2));
netcdf.putVar(myncid, varid23, prctile(tasmax_ens_opt(21:50,:),95,2));
netcdf.close(myncid);
%**************************************************************************
clear all 
clc

lat = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lon');
sst = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','sst'); % since 1854
sst = sst(:,:,1273:2016); % 1960-2021
sst(find(sst<-1e+36)) = NaN;

sst = reshape(sst,[180*89 12 62]);
sst = sst - repmat(mean(sst,3),[1 1 62]); % anomaly

sst_ann = squeeze(mean(sst(:,4:9,:),2)); % fire season

nansum = sum(isnan(sst_ann),2);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*89 1]);
lat1d  = reshape(lat2d,[180*89 1]);
weight = cos(pi*lat/180.0);
weight = repmat(weight',[180 1]);
weight = reshape(weight,[180*89 1]);

lct_1 = find(lon1d>=140&lon1d<=215&lat1d>=25&lat1d<=45&nansum==0);
weights_1 = weight(lct_1)/mean(weight(lct_1));

lct_2 = find(lon1d>=170&lon1d<=270&lat1d>=-10&lat1d<=10&nansum==0);
weights_2 = weight(lct_2)/mean(weight(lct_2));

lct_3 = find(lon1d>=150&lon1d<=200&lat1d>=-50&lat1d<=-15&nansum==0);
weights_3 = weight(lct_3)/mean(weight(lct_3));

sst_ann_1 = mean(sst_ann(lct_1,:).*repmat(weights_1,[1 62]),1);
sst_ann_2 = mean(sst_ann(lct_2,:).*repmat(weights_2,[1 62]),1);
sst_ann_3 = mean(sst_ann(lct_3,:).*repmat(weights_3,[1 62]),1);

ipo = sst_ann_2 - 0.5*(sst_ann_1+sst_ann_3);
ipo  =  movmean(ipo,9);

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');

pre = pre(:,:,59*12+1:1452); % 1960-2021
pre = reshape(pre,[720*360 12 62]);
pre = pre- repmat(mean(pre,3),[1 1 62]); % anomaly
pre = reshape(pre,[720*360 62*12]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);
pre_avg = mean(pre.*repmat(weights,[1 62*12]),1);
pre_avg  = reshape(pre_avg,[12 62]);
pre_ann = sum(pre_avg(4:9,:),1);

[R,P] = corrcoef(ipo,pre_ann)

myncid  = netcdf.create('figs9.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'yr', 62);
varid1  = netcdf.defVar(myncid, 'ipo', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_ann', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pre_ann_movavg', 'double', [dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, ipo);
netcdf.putVar(myncid, varid2, pre_ann);
netcdf.putVar(myncid, varid3, movmean(pre_ann,9));
netcdf.close(myncid);
%**************************************************************************
clear all 
clc

lat = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lat');
lon = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','lon');
sst = ncread('H:\WORKS\29-Siberian_Fire\figure3\sst.mnmean.nc','sst'); % since 1854
sst = sst(:,:,565:2016); % 1901-2021
sst(find(sst<-1e+36)) = NaN;

sst = reshape(sst,[180*89 12 121]);
sst = sst - repmat(mean(sst,3),[1 1 121]); % anomaly

sst_ann = squeeze(mean(sst(:,4:9,:),2)); % fire season

nansum = sum(isnan(sst_ann),2);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';
lon1d  = reshape(lon2d,[180*89 1]);
lat1d  = reshape(lat2d,[180*89 1]);
weight = cos(pi*lat/180.0);
weight = repmat(weight',[180 1]);
weight = reshape(weight,[180*89 1]);

lct_spna = find(lon1d>=310&lon1d<=350&lat1d>=45&lat1d<=65&nansum==0);
weights_spna = weight(lct_spna)/mean(weight(lct_spna));

lct_glb = find(lat1d>=-65&lat1d<=65&nansum==0);
weights_glb = weight(lct_glb)/mean(weight(lct_glb));
sst_ann_glb = mean(sst_ann(lct_glb,:).*repmat(weights_glb,[1 121]),1);

for i = 1:size(sst_ann,1)  
    if nansum(i) == 0
      [b,bint,r,rint,stats] = regress(sst_ann(i,:)',[ones(size(sst_ann_glb')) sst_ann_glb']);
      sst_ann_res(i,:) = sst_ann(i,:)-(b(2)*sst_ann_glb+b(1));
    else
      sst_ann_res(i,:) = NaN;
    end
end

sst_ann_res_spna = mean(sst_ann_res(lct_spna,:).*repmat(weights_spna,[1 121]),1);
sst_ann_spna = mean(sst_ann(lct_spna,:).*repmat(weights_spna,[1 121]),1);
spna  =  movmean(sst_ann_res_spna,9);

diri = 'H:\WORKS\29-Siberian_Fire\ext_fig1\';
lon = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lon');
lat = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'lat');
pre = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'pre');
stn = ncread(strcat(diri,'cru_ts4.07.1901.2022.pre.dat.nc'),'stn');

pre = pre(:,:,1:1452); % 1901-2021
pre = reshape(pre,[720*360 12 121]);
pre = pre- repmat(mean(pre(:,:,60:end),3),[1 1 121]); % anomaly
pre = reshape(pre,[720*360 121*12]);

stn = stn(:,:,1:1452); % 1960-2021
stn = reshape(stn,[720*360 12*121]);

weight = cos(pi*lat/180.0);
weight = repmat(weight',[720 1]);
weight = reshape(weight,[720*360 1]);

lon2d  = repmat(lon,[1 length(lat)]);
lat2d  = repmat(lat,[1 length(lon)])';

lon1d  = reshape(lon2d,[720*360 1]);
lat1d  = reshape(lat2d,[720*360 1]);

lct = find(lat1d>=55&lat1d<=67.5&lon1d>=95&lon1d<=125);
weights = weight(lct)/mean(weight(lct));

pre = pre(lct,:);
stn = stn(lct,:);

for i = 1:size(stn,2)
    lct2 = find(stn(:,i)>=1);
    weights2 = weights(lct2)/mean(weights(lct2));
    pre_avg(i) = mean(pre(lct2,i).*weights2,1);
    clear lct2 weights2
end

pre_avg  = reshape(pre_avg,[12 121]);
pre_ann = sum(pre_avg(4:9,:),1);

[R,P] = corrcoef(spna,pre_ann)

myncid  = netcdf.create('figs11.nc', 'NC_NOCLOBBER');
dimid1  = netcdf.defDim(myncid, 'yr', 121);
varid1  = netcdf.defVar(myncid, 'spna', 'double', [dimid1]);
varid2  = netcdf.defVar(myncid, 'pre_ann', 'double', [dimid1]);
varid3  = netcdf.defVar(myncid, 'pre_ann_movavg', 'double', [dimid1]);
netcdf.endDef(myncid);
netcdf.putVar(myncid, varid1, spna);
netcdf.putVar(myncid, varid2, pre_ann);
netcdf.putVar(myncid, varid3, movmean(pre_ann,9));
netcdf.close(myncid);






