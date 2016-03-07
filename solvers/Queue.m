classdef Queue < handle
    % 
    % FIFO (first in - first out) Queue implementation based on singly
    % linked list 
    % 
    %
    %
    

    properties (Access = private)
        
        back        % pointer to the back Node of the queue. 
        front       % pointer to the front Node of the queue.
        size        % store a size variable to keep track on how many 
                    % elements are currently queued
    end
    
    
    methods
        
        %when you create the quees
        function obj = Queue()
            
            obj.back = -1;
            obj.front = -1; 
            obj.size = 0;
    
        end
        
        
        function s = getSize(obj) 
            s = obj.size;
        end
        
        function b = isempty(obj)   % return true when the queue is empty
            b = ~logical(obj.size());
        end
        
        function s = empty(obj) % clear all the data in the queue and return the number of elements cleared! 
            s = obj.size();
            %lose the pointers to all elements...
            obj.back = -1;
            obj.front = -1;
            obj.size = 0;
        end
        
        function enqueue(obj, data)
            %package the data in a node class
            node = Node(data);
            
            if(obj.back == -1)
                obj.back = node;
                obj.front = node;
            else
                old_back = obj.back;
                
                obj.back = node;
                node.setNext(old_back);
                old_back.setPrev(node);
            end
            obj.size = obj.size + 1;
        end
        
        function b = getBackNode(obj)
            
            b = obj.back;
        
        end
        
        function f = getFrontNode(obj)
            f = obj.front;
        end
       
        
        function data = dequeue(obj) % Ԫ
            
            if obj.size == 0 
                error('Queue:NO_Data', 'Trying to dequeue an empty queue!!');
            else
                    % if there is only 1 element left in the queue, return
                    % it and empty the queue.
                
                if(obj.size == 1) 
                    node = obj.front;
                    data = node.getData();
                    obj.empty();
                else % do a proper removal... 
                    node = obj.front;
                    data = node.getData();
                    obj.front = node.getPrev();
                    obj.front.setNext(-1);
                    obj.size = obj.size - 1;
                end
            end     
            
            
        end
        
        
        function displayForward(obj) 

            if obj.size()
                
                node = obj.back;
                k = 1;
                while( node ~= -1)
                    disp([num2str(k) '-th element of the stack:']);
                    disp(node.getData);
                    node = node.getNext();
                    k = k+1;
                end
                
            else
                disp('The queue is empty');
            end
        end
        
        function displayBackward(obj) 

            if obj.size()
                node = obj.front;
                k = 1;
                while( node ~= -1)
                    disp([num2str(k) '-th element of the stack:']);
                    disp(node.getData);
                    node = node.getPrev();
                    k = k+1;
                end
            else
                disp('The queue is empty');
            end
        end
        
 
    end
end