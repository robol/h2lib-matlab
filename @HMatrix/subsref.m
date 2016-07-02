function b=subsrefs(a,s)
  if s.type ~= '()'
    error('Unsupported selection method for HMatrix');
  end

  indices = s.subs;
  if length(indices) ~= 2
    error('HMatrices have exactly two dimensions');
  end

  ii = indices{1};
  jj = indices{2};

  n = matrix_size(a);
  n = n(2);

  b = zeros(length(ii), length(jj));

  for j = 1 : length(jj)
    ej = [ zeros(jj(j)-1,1) ; 1 ; zeros(n - jj(j), 1) ];
    x = a * ej;
    b(:,j) = x(ii);
  end

end

