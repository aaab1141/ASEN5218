function [] = save_fig_png(filename)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% This function saves a figure both at a matlab .fig file and a .png file
% so that you don't have to do it manually cause that's annoying when you
% have a million figures. The .fig files end up in a folder in the current
% folder called "figs" and the .png files end up in a folder in the main
% folder called "pngs".
% 
% INPUT: string with the filename you want your files to have
% 
% Written 2019-11-14 | Aaron Aboaf
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

curfol = cd;

if ~exist([curfol,'\Figs'], 'dir')
    mkdir([curfol,'\Figs'])
end
if ~exist([curfol,'\PNGs'], 'dir')
    mkdir([curfol,'\PNGs'])
end

cd([curfol,'\Figs'])
saveas(gcf,[filename,'.fig'])

cd([curfol,'\PNGs'])
saveas(gcf,[filename,'.png'])

cd(curfol)
end

