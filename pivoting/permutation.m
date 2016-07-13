function P = permutation(vec_star)
% Returns a permutation matrix given a unique vector!
vec = [1 : length(vec_star)]';
P = bsxfun(@eq, vec', vec_star); P = P';
end