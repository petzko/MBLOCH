classdef NodeGPU < handle
    %NODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private) 
        content;
        next;
        prev;
    end
    
    methods
        function obj = NodeGPU(data)
            %put the data on the GPU
            obj.content = gpuArray(data); 
            obj.next = -1;
            obj.prev = -1; 
        end
        
        %return the node pointed to by next
        function n = getNext(obj)
            n = obj.next;
        end
        % set the node pointed to by next
        function n = setNext(obj,node)
            n = obj.next;
            obj.next = node;
        end
        
        
        %return the node pointed to by prev
        function n = getPrev(obj)
            n = obj.prev;
        end
        % set the node pointed to by prev
        function n = setPrev(obj,node)
            n = obj.prev;
            obj.prev = node;
        end
        
        function data = getData(obj)
           data = obj.content; 
        end
        
    end
    
end

