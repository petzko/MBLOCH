classdef Controller < handle
    %CONTROLLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        viewhandle; 
        sim_settings; 
        job; 
        x;
    
    end
    
   
    methods ( Access = public )
        
        function obj = Controller(view)
            display('controller initialized');
            obj.viewhandle = view; 
            
            obj.x = [1,2,3];
        end
        
        function filename = openfile(obj,srcHandle,eventData)
            display('open file launched');
            filename = uigetfile('*.set');
            obj.sim_settings = parseInput(filename);
        end

        function newsim(obj,srcHandle,eventData)
            display('setup new sim menu launched');
        end
        
        function saveworkspacedata(obj,srcHandle,eventData)
            display('save data launched');
            display(obj.x);
            x = obj.x;
            save('testsave.mat','x');
            
         end
        
        function startsim(obj,srcHandle,eventData)
             
            if(obj.view_handle.init > 0);
                display('init > 0');
            else
                display('init < 0')
            end
            
        end
  
    end
    
    
end

