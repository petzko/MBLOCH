function [ output_args ] = dlegend( legInf,legTitle )


    hLegend =legend(legInf);
    set(hLegend,'FontName','Arial','FontSize',10.0,'Location','northeast'); 
    hlt = text(...
        'Parent', hLegend.DecorationContainer, ...
        'String', legTitle, ...
        'FontName','Arial',...
        'FontSize',12.0,...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Position', [0.5, 1.05, 0], ...
        'Units', 'normalized');

end

