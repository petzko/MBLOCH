function drawnowexample(option)

persistent f
if ~nargin

f = figure('doublebuffer','on');
b = uicontrol('parent',f,'callback','drawnowexample(''kill'')','string','kill!')


while 1
    p = rand(1,1000000);
    pl = plot(p(1:100:1000000));
    drawnow
    if ~ishandle(f)
        break
    end
end

elseif strcmp(option,'kill')
    
    close(f)
end
