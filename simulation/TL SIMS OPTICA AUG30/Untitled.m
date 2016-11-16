% fh = gcf; 
% nplots = length(fg.Children);
% xdat = cell(nplots,1);
% ydat = cell(nplots,1);
% for p = 1:nplots
%     dat = fh.Children(p);
%     xdat_ = dat.Children.XData;
%     idx = (xdat_ >= dat.XLim(1) & xdat_ <= dat.XLim(2));
%     ydat{p} = dat.Children.YData(idx);
%     xdat{p} = dat.Children.XData(idx);
% end
% 
% % test... 
nplots = 6
for p = 1:nplots
    subplot(nplots,1,p); 
    plot(xdat{p},ydat{p})
end