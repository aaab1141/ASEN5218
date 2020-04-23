function [] = save_fig_png(filename)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This function saves a figure both at a matlab .fig file and a .png file
% so that you don't have to do it manually cause that's annoying when you
% have a million figures. The .fig files end up in a folder in the Plots
% folder called "Figs", the .png files end up in a folder in the Plots
% folder called "PNGs", and the .pdf files end up in a folder in the Plots 
% folder called "PDFs".
% 
% INPUT: string with the filename you want your files to have
% 
% Written 2019-11-14 | Aaron Aboaf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

curfol = cd;

if ~exist([curfol,'\Plots'],'dir')
    mkdir([curfol,'\Plots'])
end
if ~exist([curfol,'\Plots\Figs'],'dir')
    mkdir([curfol,'\Plots\Figs'])
end
if ~exist([curfol,'\Plots\PNGs'], 'dir')
    mkdir([curfol,'\Plots\PNGs'])
end
if ~exist([curfol,'\Plots\PDFs'], 'dir')
    mkdir([curfol,'\Plots\PDFs'])
end

cd([curfol,'\Plots\Figs'])
saveas(gcf,[filename,'.fig'])

cd([curfol,'\Plots\PNGs'])
print(gcf,[filename,'.png'],'-dpng','-r350')

cd([curfol,'\Plots\PDFs'])
saveas(gcf,[filename,'.pdf'])


cd(curfol)
end

