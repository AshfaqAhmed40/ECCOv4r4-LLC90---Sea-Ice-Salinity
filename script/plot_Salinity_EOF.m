clear all; close all; clc;

%% Data
% tic % 200 seconds
% HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
% cd(strcat(HOME,'input/data/'));
% sourceDir = (strcat(HOME,'input/data/'));
% files = dir(fullfile(sourceDir, '*.nc'));
% SALT_3D = zeros(90, 90, length(files)); % preparing a 3D variable

% 
% for ii = 1:length(files)
%     fname = files(ii).name;
%     readncfile;
%         clear A i i_g j j_g k k_l k_p1 k_u ll XC_bnds ...
%         THETA XC YC nYG YC_bnds Z_bnds Zl Zp1 Zu XG YG
%     tile = 7; % arctic slice
%     depth = 1; % upper layer only
%     SALT = squeeze(SALT(:,:,tile,depth)); % considering Arctic slice
%     SALT = flip(SALT);
%     SALT_3D(:,:,ii) = SALT;
% end
% toc

HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input/data/'));
load SALT_3D.mat
load coords.mat
 clc;
 
% Create a datetime vector for each month from 1992 to 2017
startDate = datetime(1992, 1, 1);  % Start in January 1992
endDate = datetime(2017, 12, 1);   % End in December 2017
months = dateshift(startDate, 'start', 'month', 0:311);
obsv = datetime(months);

%% Regrid all the variables (3D to 2D)
SALTr = reshape(SALT_3D,length(SALT_3D(:,1,1))*length(SALT_3D(1,:,1)),length(SALT_3D(1,1,:)))';
latr = reshape(LAT,1,length(LAT(:,1))*length(LAT(1,:)));
lonr = reshape(LON,1,length(LON(:,1))*length(LON(1,:)));
SALTr = [latr;lonr;SALTr]; 

%% 3. Remove nan data (actually, 100% NaN columns)
SALTr(:,sum(isnan(SALTr(3:end,:)),1)==length(SALTr(3:end,1)))=[];

%% 4. Take out LAT LON
coords = SALTr(1:2,:);
SALTr(1:2,:)   =[];

%% 5. Fill in with the nearest temporal data 
SALTr = fillmissing(SALTr,'linear',1,'EndValues','near');
SSS = SALTr;
%SSS = rmmissing(SSS,2);

%% 6. Normalize the temp data
%datms = SSS./nanstd(SSS(:));
datms = SSS;

%% Singular Value Decomposition
[U,S,V] = svds([datms],5);

%% Calculate EEOFs
eeof1s = V(1:length(datms(1,:)),1);
eeof2s = V(1:length(datms(1,:)),2);
eeof3s = V(1:length(datms(1,:)),3);

%% Reshape EEOFs
finaleeof1s = nan(size(LAT,1),size(LAT,2));
finaleeof2s = nan(size(LAT,1),size(LAT,2));
finaleeof3s = nan(size(LAT,1),size(LAT,2));

for ff = 1:length(eeof1s(:,1))
    finaleeof1s(coords(1,ff)==LAT(:,1),coords(2,ff)==LON(1,:)) = eeof1s(ff,1);
    finaleeof2s(coords(1,ff)==LAT(:,1),coords(2,ff)==LON(1,:)) = eeof2s(ff,1);
    finaleeof3s(coords(1,ff)==LAT(:,1),coords(2,ff)==LON(1,:)) = eeof3s(ff,1); 
end

%% 10. Saving the variables
U1 = U(:,1); U2 = U(:,2); U3 = U(:,3);

%% 11. Prepare the color map
clr = flipud(brewermap(15,'BrBG'));
clri= interp1(1:1:15,clr,1:0.25:15,'linear');


%% Plot all three EEOFs
clc;
EOF1 = (diag(S(1,1)))^2/(sum(diag(S)))^2
EOF2 = (diag(S(2,2))+diag(S(1,1)))^2/(sum(diag(S)))^2
EOF3 = (diag(S(2,2))+diag(S(1,1))+diag(S(3,3)))^2/(sum(diag(S)))^2

figure(1),
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Plot the EOF1
nexttile;
[~, hContour]=contourf(rot90(rot90(finaleeof1s).*U1(1,1)),200,'edgecolor','none');
set(gca,'Color',[0.45 0.45 0.45]); shading interp;
colormap(gca,clri); caxis([-1.5e-3 1.5e-3]);
set(gca, 'XTickLabel',[]); set(gca, 'YTickLabel',[]);
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');

      
% Plot the EOF2
nexttile
[~, hContour]=contourf(rot90(rot90(finaleeof2s).*U2(1,1)),200,'edgecolor','none');
set(gca,'Color',[0.45 0.45 0.45]); shading interp;
colormap(gca,clri); caxis([-1.5e-3 1.5e-3]);
set(gca, 'XTickLabel',[]); set(gca, 'YTickLabel',[]);
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');
        
% Plot the EOF3
nexttile
[~, hContour]=contourf(rot90(rot90(finaleeof3s).*U3(1,1)),200,'edgecolor','none');
set(gca,'Color',[0.45 0.45 0.45]); shading interp;
colormap(gca,clri); caxis([-1.5e-3 1.5e-3]);
set(gca, 'XTickLabel',[]); set(gca, 'YTickLabel',[]);
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3);
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');
    
c = colorbar;
c.FontSize = 25;
c.FontWeight = 'bold'; c.Label.String = 'Loadings';
c.Layout.Tile = 'south';
set(gcf, 'Position', [-2550,616,1208,451]);


%% Test for nice picture
% figure(100),
% clf;
% t = tiledlayout(6,5);
% t.TileSpacing = 'compact';
% for i = 31:60
% nexttile
% [~, hContour]=contourf(rot90(rot90(finaleeof3s).*U3(i,1)),200,'edgecolor','none');
% shading interp; colormap(gca,clri1);
% caxis([-1.5e-3 1.5e-3]);
% title (sprintf('%d',i));
% set(gca, 'YTickLabel',[]); set(gca, 'XTickLabel',[]);
% set(gcf, 'Position', [-2135,61,914,1276]);
% end

 %% Animation of the EOFs
% clc; close all;
% choice = 3;
% switch choice
%     case 1
%         finaleeof = finaleeof1s;
%         U = U1;
%         eofFolderName = 'EOF1-Salinity';
%     case 2
%         finaleeof = finaleeof2s;
%         U = U2;
%         eofFolderName = 'EOF2-Salinity';
%     case 3
%         finaleeof = finaleeof3s;
%         U = U3;
%         eofFolderName = 'EOF3-Salinity';
%     otherwise
%         error('Invalid choice');
% end
% 
% outputFolder = fullfile(HOME, 'figures', 'EOF',  'Animation', eofFolderName);
% cd(outputFolder)
% 
% for i=1:100
%     fig = figure('Visible', 'off');
%     [~, hContour] = contourf(rot90(rot90(finaleeof).*U(i,1)), 25, 'edgecolor', 'k');
%     set(gca,'Color',[0.45 0.45 0.45]); shading interp;
%     colormap(gca,clri); caxis([-1.5e-3 1.5e-3]);
%     set(gca, 'XTickLabel',[]); set(gca, 'YTickLabel',[]);
%     hold on; box on;
%     set(gca, 'GridLineStyle', '-', 'LineWidth', 3);
%     set(gca, 'FontSize', 12, 'FontWeight', 'bold');
%     set(gca, 'Layer', 'top');
%     title(sprintf('%s', string(obsv(i))),'fontsize',25)
%     fileName = sprintf('%s - %d.jpg', eofFolderName, i);
%     fullFilePath = fullfile(outputFolder, fileName);
%     saveas(gcf, fullFilePath);    
% end


%% 13. Plot EEOFs
figure(2), clf;
t = tiledlayout(3,1);
t.TileSpacing = 'tight';
t.Padding = 'compact';


nexttile
plot(months,U1,'k-','markersize',1,'linewidth',3); hold on;
legend('EEOF1', 'location','northwest');
ylim([0.054 0.06]);
title('Salinity EOFs')
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3);
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');
set(gca, 'YTickLabel'); set(gca, 'XTickLabel',[]);

nexttile
plot(months,U2,'b-','markersize',1,'linewidth',3); hold on;
legend('EEOF2', 'location','northwest');
ylim([-.2 .2]);
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3);
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');
set(gca, 'YTickLabel'); set(gca, 'XTickLabel',[]);

nexttile
plot(months,U3,'r-','markersize',1,'linewidth',3); hold on;
legend('EEOF3', 'location','northwest');
ylim([-.2 .2]);
hold on; box on; 
set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3);
set(gca, 'FontSize', 20, 'FontWeight', 'bold');
set(gca, 'Layer', 'top');
set(gca, 'YTickLabel');
set(gcf, 'Position', [-2550,616,1208,451]);


%% 14. Plotting the frequency spectra
figure(3), clf;
yg(1,:) = U1;
yg(2,:) = U2;
yg(3,:) = U3;


[pxx1,f] = plomb(yg(1,:),obsv,[],10,'power'); 
[pxx2,f] = plomb(yg(2,:),obsv,[],10,'power'); 
[pxx3,f] = plomb(yg(3,:),obsv,[],10,'power'); 
f = f*86400*365;


loglog(f,10000.*pxx1,'k-','linewidth',2); hold on
loglog(f,pxx2,'b-','linewidth',2); hold on
loglog(f,pxx3,'r-','linewidth',2); hold on
set(gcf, 'Position', [-4691,599,766,593]); 
set(gca, 'GridLineStyle', '-', 'LineWidth', 3); hold on
set(gca,'XTickLabel');set(gca,'YTickLabel');
set(gca,'XTick');set(gca,'YTick');
ylabel('Magnitude'); 
xlabel('Frequency (Year^{-1})');
box on;
set(gca, 'FontSize', 16, 'FontWeight', 'bold');
xlim([10e-2 5]);
ylim([10e-16 10e4]); 
legend('EOF1','EOF2','EOF3','fontsize',16);
title('Salinity power spectrum')
