clear all; close all; clc;

%% Data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input'));
fname = "ExampleFile-Salinity-Temperature-2017-12-01.nc";
readncfile;
tile = 7; % arctic slice
depth = 1:50; Z = Z(depth); % upto 2 km
SALT = squeeze(SALT(:,:,tile,depth)); % considering Arctic slice upto 2 km
XC = squeeze(XC(:,:,tile)); YC = squeeze(YC(:,:,tile));
SALT = flip(SALT);
n = 1;

clear A fname i i_g j j_g k k_l k_p1 k_u ll ...
    XC_bnds YG YC_bnds Z_bnds Zl Zp1 Zu XG


%% Salinity
clr = flipud(brewermap(15,'BrBG'));
clri= interp1(1:1:15,clr,1:0.25:15,'linear');

figure(1),clf;
t = tiledlayout(1,3);
t.TileSpacing = 'tight';
t.Padding = 'compact';

xIdx = 40; %intersect on x-axis, parallel to y-axis
yIdx = 79; %intersect on y-axis, parallel to x-axis
intrsctd_SSS = cell(50,90,2);
intrsctd_SSS{1} = fliplr(squeeze(SALT(:,90-xIdx,:))'); %against xIdx
intrsctd_SSS{2} = fliplr(squeeze(SALT(90-yIdx,:,:))'); %against yIdx

% Plotting
nexttile
contourf(rot90(rot90(SALT(:,:,1))),100,'edgecolor','none'); shading interp; 
colormap(gca, clri); set(gca,'color',[0.55 0.55 0.55]);
line([xIdx, xIdx], [1, 90], 'Color', 'b', 'LineWidth', 4); hold on;
line([1, 90], [yIdx, yIdx], 'Color', 'r', 'LineWidth', 4); hold on;
box on; set(gca, 'GridLineStyle', '-','LineWidth', 3); hold on; 
set(gca, 'Layer', 'top');
title('Intersection line on Arctic')
xlabel('x-index'); ylabel('y-index'); set(gca,'FontSize', 15)



nexttile
contourf(intrsctd_SSS{1},100,'edgecolor','none');
colormap(gca, clri); set(gca,'color',[0.55 0.55 0.55])
set(gca, 'YDir', 'reverse');
box on; set(gca, 'GridLineStyle', '-','LineWidth', 3);
hold on; set(gca, 'Layer', 'top');
title('Salinity with depth (blue intersect)'); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel'); 
ylabel('Depth (m)'); set(gca,'FontSize', 15)

nexttile
contourf(intrsctd_SSS{2},100,'edgecolor','none');
colormap(gca, clri); set(gca,'color',[0.55 0.55 0.55])
set(gca, 'YDir', 'reverse');
box on; set(gca, 'GridLineStyle', '-','LineWidth', 3);
hold on; set(gca, 'Layer', 'top');
title('Salinity with depth (red intersect)'); 
set(gca,'XTickLabel',[]); set(gca,'YTickLabel'); 
ylabel('Depth (m)'); set(gca,'FontSize', 15)
set(gcf, 'Position', [-2348,775,1403,476]);




c = colorbar;
c.FontSize = 15; c.FontWeight = 'bold';
c.Layout.Tile = 'east'; caxis([29.5 36.5]); clc;