function varargout = dfigure(varargin)

% Optional argument DName gives figure a name and replaces the current figure
dnlocs = find(cellfun(@(x)isequal(x,'DName'),varargin));
if ~isempty(dnlocs)
    myname = varargin{dnlocs(1)+1};
    redvarargin = [varargin(1:dnlocs(1)-1),varargin(dnlocs(1)+2:end)];
    fn=findobj('Type','figure','-and','Name',myname);
    if isempty(fn)
        fn=figure('Name',myname,redvarargin{:});
    else
        fn=figure(fn(1)); if ~isempty(redvarargin), set(fn,redvarargin{:}); end
    end
else
    fn = figure(varargin{:});
end



my_FontName = 'Arial';
my_AxesLineWidth = 1.0;
my_LineWidth = 1.0;
my_Color = 'w';
my_ColorOrder =   [  0.3482    0.7424    0.5473
    1.0000         0         0
     0.5556    0.3472    0.2211
    0    0.5000         0  
    0.9880    0.8066    0.1794
    0.0972    0.0972    0.1806
    0.0723    0.4887    0.8467
    0.2081    0.1663    0.5292];

set(gcf,'Color',my_Color);
set(gcf,'DefaultAxesColorOrder',my_ColorOrder);
set(gcf,'DefaultLineLineWidth',my_LineWidth);
set(gcf,'DefaultAxesLineWidth',my_AxesLineWidth);
set(gcf,'DefaultAxesFontName',my_FontName);
set(gcf,'DefaultTextFontName',my_FontName);
set(gcf,'DefaultAxesFontSize',13.0)

% movegui(gcf,'center');

varargout = {fn};

end

