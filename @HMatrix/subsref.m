function b=subsrefs(a,s)
  indices = s.subs;
  if (max(size(a)) > 1) && (strcmp(s.type, '()'))
    b = a(indices{:});
    return;
  end

  if s.type == '.'
    if s.subs == 'row_cluster'
      b = a.row_cluster;
    end
    if s.subs == 'col_cluster'
      b = a.col_cluster;
    end

    return;
  end

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

