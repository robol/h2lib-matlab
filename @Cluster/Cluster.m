classdef Cluster < handle
    %CLUSTER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % The pointer to the C structure holding the cluster.
        cluster
    end
    
    methods
        
        function obj = Cluster(n, k)
            std_subdivision_scheme (obj, n, k);
        end
        
        function delete(obj)
            release(obj);
        end       
        
        function disp(obj)            
            fprintf('    Cluster of size %d\n\n', cluster_size(obj));
        end
    end
    
    methods (Access = private)
        function release(c)
            delete_cluster(c);
            obj.cluster = 0;
        end
    end
    
end

