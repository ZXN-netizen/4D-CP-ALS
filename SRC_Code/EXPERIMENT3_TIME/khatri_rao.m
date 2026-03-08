function KR = khatri_rao(U, V)
    R = size(U, 2); KR = zeros(size(U, 1) * size(V, 1), R);
    for r = 1:R, KR(:, r) = kron(U(:, r), V(:, r)); end
end