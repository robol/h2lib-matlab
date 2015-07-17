classdef H2Matrix < handle
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
        h2matrix
    end
    
    methods
        function obj = H2Matrix(varargin)            
            if (length (varargin) < 1)
                obj.hmatrix = 0;
                obj.row_cluster = 0;
                obj.col_cluster = 0;
            else
                kind = varargin{1};
                
                switch (kind)
                    case 'tridiagonal'
                        if (length(varargin) < 6)
                            fprintf ('You need to specify exactly 3 vectors for tridiagonal shape');
                            return;
                        end
                        tridiag (obj, varargin{2}, varargin{3}, varargin{4}, varargin{5}, varargin{6});
                        obj.row_cluster = varargin{2};
                        obj.col_cluster = varargin{3};
                    case 'band'
                        create_band_hmatrix(obj, varargin{:});
                    case 'generators'
                        if (length(varargin) < 10)
                            fprintf ('You need to specify exactly 5 vectors and 2 integers for generators shape');
                            return;
                        end
	                generators (obj, varargin{2}, varargin{3}, varargin{4}, varargin{5}', varargin{6}', varargin{7}', varargin{8}', varargin{9}, varargin{10});
                        obj.row_cluster = varargin{2};
                        obj.col_cluster = varargin{3};
                    case 'lowrank'
                        U = varargin{4};
                        V = varargin{5};
                        n = size(U,1);
                        k = size(U,2);
                        generators (obj, varargin{2}, varargin{3}, arrayfun(@(i) U(i,:)*V(i,:)', 1 : n), U', V', U', V', k, k)
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
        function create_band_hmatrix(obj, varargin)
            if (length(varargin) < 6)
                fprintf ('You need to specify exactly 3 vectors and 2 integers for band shape');
                return;
            end

            if length(varargin) == 6
                p = varargin{5};
                q = varargin{6};
                A = varargin{4};
                n = size(A,1);
                L = zeros(p,n-1);
                U = zeros(q,n-1);
                for i = 1 : p
                    L(i,1:n-i) = diag(A,-i);
                end
                for i = 1 : q
                    U(i,1:n-i) = diag(A,i);
                end
                d = diag(A);
            else
                L = varargin{5};
                U = varargin{6};
                p = varargin{7};
                q = varargin{8};
                d = varargin{4};
            end
			
            band (obj, varargin{2}, varargin{3}, d, L, U, p, q);
            obj.row_cluster = varargin{2};
            obj.col_cluster = varargin{3};            
        end

        function release(obj)
            if obj.hmatrix ~= 0
                delete_hmatrix(obj);
                obj.hmatrix = 0;
            end
        end
    end
    
end

