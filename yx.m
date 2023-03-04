function lmd=yx(ZZ,Z)

normA = sqrt(sum(ZZ .^ 2, 2));
normB = sqrt(sum(Z .^ 2, 2));
lmd = bsxfun(@rdivide, bsxfun(@rdivide, A * B, normA), normB);
end
