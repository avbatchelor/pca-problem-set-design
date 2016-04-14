function my_save_fig(filename)

% Make the figure full screen 
p = get(0,'MonitorPositions');
set(gcf,'Position',p(1,:)) % Makes the figure full screen on the first monitor

% Set figure size 
if ~exist('figSize','var')
    figSize = [5 3];
end
set(gcf, 'PaperSize', figSize);
set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0,figSize]);

% Save as emf file 
imageFilename = [filename,'.emf'];
print(gcf,'-dmeta',imageFilename,'-r300','-painters','-cmyk')
