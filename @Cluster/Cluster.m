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
        
        function release(c)
            delete_cluster(c);
        end
    end
    
end

