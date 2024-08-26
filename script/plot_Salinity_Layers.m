clear all; close all; clc;

%% Data
HOME = '/Users/aahmed78/Desktop/PhD/Arctic Ocean Salinity/Ocean Model/ECCOv4r4 LLC0090/';
cd(strcat(HOME,'input'));
fname = "ExampleFile-Salinity-Temperature-2017-12-01.nc";
readncfile;
tile = 7; % arctic slice

Z = round(Z,2);
depth = 1:37; Z = Z(depth); % upto 2 km
SALT = squeeze(SALT(:,:,tile,depth)); %taking only the Arctic slice upto 2 km
SALT = flip(SALT);

clear A fname i i_g j j_g k k_l k_p1 k_u ll ...
    XC_bnds YG YC_bnds Z_bnds Zl Zp1 Zu XG


%% Thgermohaline profile
% Fram Strait (80,53)
FramStrait_SSS = squeeze(SALT(80,53,:));

% Barents Sea (53,27)
BarentsSea_SSS = squeeze(SALT(53,27,:));

% Beaufort Gyre (3,66)
BeaufortGyre_SSS = squeeze(SALT(3,66,:));

%% Salinity Layers
clr = flipud(brewermap(15,'BrBG'));
clri= interp1(1:1:15,clr,1:0.25:15,'linear');


figure(1),clf;
n = 1;
t = tiledlayout(sqrt(n),sqrt(n));
t.TileSpacing = 'tight';
t.Padding = 'compact';

for i = 1:n
   nexttile
   contourf(rot90(rot90(SALT(:,:,i))),150,'edgecolor','none'); shading interp; 
   colormap(gca, clri); caxis([31 35]);
   set(gca,'color',[0.55 0.55 0.55])
   title(sprintf('Depth = %0.2f m',Z(i)),'FontSize', 25)
   set(gcf, 'Position', [-2333,75,1372,1262]);
   set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':', 'LineWidth', 3);
   box on; set(gca, 'GridLineStyle', '-','LineWidth', 3);
   hold on; set(gca, 'Layer', 'top');
   set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
   
   % plot the grid
        cellsToHighlight = [80, 53; 53, 27; 3, 66];
        hold on;
        for idx = 1:size(cellsToHighlight, 1)
            newRow = size(SALT,1) - cellsToHighlight(idx, 1) + 1;
            newCol = size(SALT,2) - cellsToHighlight(idx, 2) + 1;
            xCoords = [newCol-1.5, newCol+1.5, newCol+1.5, newCol-1.5];
            yCoords = [newRow-1.5, newRow-1.5, newRow+1.5, newRow+1.5];
            patch(xCoords, yCoords, 'y'); 
        end
        hold off;


end

c = colorbar;
c.FontSize = 30;
c.FontWeight = 'bold';
c.Layout.Tile = 'east';