classdef HMatrix < handle
    %HMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Row cluster object. This is held here to make sure
        % that we keep a reference of it around, preventing its
        % deletion. 
        row_cluster
        
        % Column cluster object. The same of above holds. 
        col_cluster
        
        % Pointer to the hmatrix struct in C
        hmatrix
    end
    
    methods
        function obj = HMatrix(varargin)            
            if (length (varargin) < 1)
                obj.hmatrix = 0;
                obj.row_cluster = 0;
                obj.col_cluster = 0;
            else
                kind = varargin{1};
                
                switch (kind)
                    case 'tridiagonal'
                        if (length(varargin) < 4)
                            fprintf ('You need to specify exactly 3 vectors for tridiagonal shape');
                            return;
                        end
                        tridiag (obj, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6});
                        obj.row_cluster = varargin{2};
                        obj.col_cluster = varargin{3};
                    case 'pointer'
                        obj.hmatrix = varargin{2};
                        obj.row_cluster = varargin{3};
                        obj.col_cluster = varargin{4};
                    otherwise
                        fprintf ('Unsupported matrix kind specified');
                end
            end
        end
        
        function delete(obj)
            release(obj);
        end
        
        function disp(obj)
            sz = matrix_size(obj);
            fprintf('    HMatrix of size %d x %d, rank %d\n\n', ...
                sz(1), sz(2), hmatrix_rank (obj));
            disp(full(obj));
        end
        
    end
    
    methods (Access = private)
        function release(obj)
            if obj.hmatrix ~= 0
                delete_hmatrix(obj);
                obj.hmatrix = 0;
            end
        end
    end
    
end

