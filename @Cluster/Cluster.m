classdef Cluster < handle
    %CLUSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The pointer to the C structure holding the cluster.
        cluster
        
        % If true it means that the indexes of this cluster 
        % will be freed as soon as it is unreferenced. 
        owned
        
        % Pointer to the parent Cluster, if any
        parent
    end
    
    methods

	function obj = Cluster(varargin)
            if (length(varargin) == 0)
                obj.cluster = 0;
                obj.owned = false;
                obj.parent = 0;
            else               
                n = varargin{1};
                k = varargin{2};
                std_subdivision_scheme (obj, n, k);
                obj.owned = true;
                obj.parent = 0;
            end
        end   
        
        function delete(obj)
            if obj.owned
                release(obj);
            else
                % Unreference the parent, so if we are the last
                % reference MATLAB can release it. 
                obj.parent = 0;
            end
        end       
        
        function disp(obj)            
            if obj.cluster ~= 0
                fprintf('    Cluster of size %d\n\n', ...
                        cluster_size(obj));
            else
                fprintf('    Empty cluster\n\n');
            end
        end
    end
    
    methods (Access = private)
        function release(c)
            delete_cluster(c);
            obj.cluster = 0;
        end
    end
    
end

