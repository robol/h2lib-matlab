classdef HMatrix < handle
    %HMATRIX Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Pointer to the hmatrix struct in C
        hmatrix
    end
    
    methods
        function obj = HMatrix(varargin)
            if (length (varargin) < 1)
                obj.hmatrix = 0;
            else
                kind = varargin{1};
                
                switch (kind)
                    case 'tridiagonal'
                        if (length(varargin) < 4)
                            fprintf ('You need to specify exactly 3 vectors for tridiagonal shape');
                            return;
                        end
                        tridiag (obj, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6});
                    case 'pointer'
                        obj.hmatrix = varargin{2};
                    otherwise
                        fprintf ('Unsupported matrix kind specified');
                end
            end
        end
        
        function release(obj)
            delete_hmatrix(obj);
            obj.hmatrix = 0;
        end
    end
    
end

