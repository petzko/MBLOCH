classdef View  < handle
    
    
    properties
        
        % ###########################################################
        %GUI items
        Figure;
        Axes;
        Statusbar; 
        
        
        
        LoopControlPanel;
        Startbtn;Stopbtn;Pausebtn;Resumebtn;
        
        SolverSelectPanel;
        bg; %button group
        rnfdradio;laxradio; msstepstext;nrstepsedit;
        
        %menu window
        filemenu;
        newsim_menuitm; open_menuitm; save_menuitm; close_menuitm;
        % ###########################################################
        %end of GUI items!!
      
        filename;
        init;
        controller;
        title; %title of the main window
        dim; % xy dimension of the main window
        x0; %x-coordinate of the lower left corner of the figure;
        y0; %y-coordinate of the lower left corner of the figure;
        margin; %in percents
    end
    methods
        function obj = View(varargin)
            assert(length(varargin) <= 2,'error. cannot create view! too many input arguments');
            
            
            obj.title = 'Reduced Maxwell-Bloch equations simulateor';
            obj.dim = [600,500];
            obj.x0 = 100 ; obj.y0 = 100;
            
            obj.margin = 0.05;
            obj.init = 1; 
            
            %handle initial user input
            for i = 1:length(varargin)
                argi = varargin{i};
                if(strcmp(argi{1},'title'))
                    obj.title = argi{2};
                else if (strcmp(argi{1},'size') )
                        obj.dim = argi{2};
                    end
                end
            end
            
            %Initialize the controller
            obj.filename = 'none'
            obj.controller = Controller(obj);
            
            %now setup the gui;
            
            %initialize the main window
            obj.Figure = figure('Visible','on','MenuBar','None',  'Position',[obj.x0,obj.y0,obj.x0+obj.dim(1),obj.y0+obj.dim(2)]);
            %initialize the menus
            obj.filemenu = uimenu(obj.Figure,'Label','File');
            obj.newsim_menuitm = uimenu(obj.filemenu,'Label','Setup new simulation file','Callback',@obj.newsim);
            obj.open_menuitm = uimenu(obj.filemenu,'Label','Open simulation file','Callback',@obj.openfile);
            obj.save_menuitm = uimenu(obj.filemenu,'Label','Save simulation data','Callback',@obj.savedata);
            obj.close_menuitm = uimenu(obj.filemenu,'Label','Close','Callback',@obj.close);
            
            %initialize the axis;
            obj.Axes = axes(...    % Axes for plotting the selected plot
                'Parent', obj.Figure, ...
                'Units', 'normalized', ...
                'HandleVisibility','callback', ...
                'Position',[obj.margin 3*obj.margin 0.60 0.70]);
            obj.Statusbar = uicontrol(obj.Figure,'Style','text','Units','normalized', ...
                'Position',[obj.margin 0 (1-2*obj.margin) obj.margin],'FontSize',9,'HorizontalAlignment','left','BackgroundColor',[71 148 148]/255);
            obj.Statusbar.String  = ['Simulation file:' obj.filename] ;
            %initialize the loop control panel and loop control buttons
            pos = obj.Axes.Position;
            pos(2) = obj.margin + pos(2)+pos(4);
            obj.LoopControlPanel = uipanel('Parent',obj.Figure,'Position', [pos(1) pos(2) pos(3) 0.1 ], 'Title','Loop control');
            
            %Startbtn;Stopbtn;Pausebtn;Resumebtn;
            
            obj.Startbtn = uicontrol(obj.LoopControlPanel,'Style','pushbutton','String','Start','Units','normalized','Position',[0.05 0.1 0.2 0.7],'Callback',@obj.startsim);
            obj.Stopbtn = uicontrol(obj.LoopControlPanel,'Style','pushbutton','String','Stop','Units','normalized','Position',[0.35 0.1 0.2 0.7 ],'Callback',@obj.stopsim);
            obj.Pausebtn = uicontrol(obj.LoopControlPanel,'Style','pushbutton','String','Pause','Units','normalized','Position',[0.50 0.1 0.2 0.7 ],'Callback',@obj.pausesim);
            obj.Resumebtn = uicontrol(obj.LoopControlPanel,'Style','pushbutton','String','Resume','Units','normalized','Position',[0.75 0.1 0.2 0.7 ],'Callback',@obj.resumesim);
            align([ obj.Startbtn  obj.Stopbtn  obj.Pausebtn  obj.Resumebtn],'distribute','bottom');
            
            %initialize the solver settings panel and solver settings  components
            pos = obj.Axes.Position;
            pos(1) = pos(1) + pos(3) + obj.margin;
            pos(2) = pos(2)+0.5*pos(4);
            pos(3) = 1. - 3*obj.margin - pos(3);
            pos(4) = 0.5*pos(4);
            obj.SolverSelectPanel = uipanel('Parent',obj.Figure,'Position', [pos(1) pos(2) pos(3) pos(4) ], 'Title','Solver selector');
            
            obj.bg = uibuttongroup(obj.SolverSelectPanel,'Visible','on','Units','normalized',...
                'Position',[0.1 0.6 0.8 0.4],...
                'SelectionChangedFcn',@obj.rbselection);
            
            % Create three radio buttons in the button group.
            obj.rnfdradio = uicontrol(obj.bg,'Style',...
                'radiobutton',...
                'String','RNFD ',...
                'Units','normalized','Position',[0.1 0.55 0.9 .35 ],...
                'HandleVisibility','off');
            obj.laxradio = uicontrol(obj.bg,'Style',...
                'radiobutton',...
                'String','Lax-Wendroff ',...
                'Units','normalized','Position',[0.1 0.1 0.9 .35 ],...
                'HandleVisibility','off');
            align([ obj.rnfdradio obj.laxradio],'Center','distribute');
            
            obj.msstepstext = uicontrol(obj.SolverSelectPanel,'Style',...
                'text',...
                'String', 'MS order:1',...
                'Units','normalized','Position',[0.1 0.3 0.9 .1 ],...
                'HandleVisibility','off');
            
            obj.nrstepsedit = uicontrol(obj.SolverSelectPanel,'Style',...
                'slider',...
                'Units','normalized','Position',[0.1 0.1 0.8 .11 ],...
                'HandleVisibility','off','Min',1,'Max',5,...
                'SliderStep',[1 1]./(5-1),'Value',1,'Callback',@obj.sliderchange);
            
            align([ obj.bg obj.msstepstext  obj.nrstepsedit],'left','distribute');
        end
        
        function savefig(obj,filename)
            savefig(obj.Figure,filename);
        end
        
    end
    methods (Access = private)
        
        function newsim(obj,srcHandle,eventData)
            obj.controller.newsim(srcHandle,eventData);
        end
        
        
        function openfile(obj,srcHandle,eventData)
           obj.filename =  obj.controller.openfile(srcHandle,eventData);
           obj.Statusbar.String  = ['Simulation file:' obj.filename] ;
        end
        
        
        function savedata(obj,srcHandle,eventData)
            obj.controller.saveworkspacedata(srcHandle,eventData);
        end
        
        
        function close(obj,srcHandle,eventData)
            close(obj.Figure); 
        end
        
        function startsim(obj,srcHandle,eventData)
            obj.init = +1;
            obj.controller.startsim(srcHandle,eventData);
            obj.switchocntrols('off');
        end
        
        function switchocntrols(obj,mode)
        
            if(strcmp(mode,'on'))
                obj.laxradio.Enable = 'on';
                obj.rnfdradio.Enable = 'on';
                obj.nrstepsedit.Enable = 'on';
            else
                obj.laxradio.Enable = 'off';
                obj.rnfdradio.Enable = 'off';
                obj.nrstepsedit.Enable = 'off';
            end
        
            
        end
        
        function stopsim(obj,srcHandle,eventData)
           obj.init = 1;
           obj.controller.stopsim(srcHandle,eventData);
           obj.switchocntrols('on');
        end
        
        function pausesim(obj,srcHandle,eventData)
            obj.init = -1;
            obj.controller.stopsim(srcHandle,eventData);
        end
        
        function resumesim(obj,srcHandle,eventData)
            obj.init = -1; 
            obj.controller.startsim(srcHandle,eventData);            
        end
        
        function rbselection(obj,srcHandle,eventData)
            
        end
        
        
        function sliderchange(obj,srcHandle,eventData)
            
            obj.msstepstext.String = ['MS order:' num2str(srcHandle.Value)];
        end
        
    end
    
    
end
