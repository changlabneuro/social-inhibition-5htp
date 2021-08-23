function [ps, labs, each_I] = permute_slope_differences(X, Y, labels, iters, each, compare)

[labs, each_I] = keepeach_or_one( labels', each );
ps = nan( numel(each_I), 1 );

parfor i = 1:numel(each_I)
  ind = each_I{i};
  
  compare_I = findall( labels, compare, ind );
  assert( numel(compare_I) == 2, 'Expected 2 groups to compare; got %d.' ...
    , numel(compare_I) );
  
  num_a = numel( compare_I{1} );
  null_diffs = zeros( iters, 1 );
  
  [lm1, lm2] = fit_models( X, Y, compare_I{1}, compare_I{2} );
  real_diff = slope_comparison( lm1, lm2 );
  
  for j = 1:iters
    subset_a = sort( ind(randperm(numel(ind), num_a)) );
    subset_b = setdiff( ind, subset_a );
    
    [lm1, lm2] = fit_models( X, Y, subset_a, subset_b );
    null_diffs(j) = slope_comparison( lm1, lm2 );
  end
  
  ps(i) = sum( null_diffs > real_diff ) / iters;
end

end

function [lm1, lm2] = fit_models(amp, vel, ind1, ind2)

a1 = amp(ind1);
v1 = vel(ind1);

a2 = amp(ind2);
v2 = vel(ind2);

lm1 = fitlm( a1, v1 );
lm2 = fitlm( a2, v2 );

end

function b = slope_comparison(lm1, lm2)

beta1 = lm1.Coefficients.Estimate(2);
beta2 = lm2.Coefficients.Estimate(2);

b = abs( beta1 - beta2 );

end