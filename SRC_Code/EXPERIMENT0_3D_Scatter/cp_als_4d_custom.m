function [A, B, C, S] = cp_als_4d_custom(Y, R, max_iter, tol, A_t, B_t, C_t, S_t)
    [I, J, M, K] = size(Y);
    Y1 = reshape(Y, I, J*M*K); Y2 = reshape(permute(Y, [2, 1, 3, 4]), J, I*M*K);
    Y3 = reshape(permute(Y, [3, 1, 2, 4]), M, I*J*K); Y4 = reshape(permute(Y, [4, 1, 2, 3]), K, I*J*M);
    A = A_t + 0.05*(randn(I, R) + 1j*randn(I, R)); B = B_t + 0.05*(randn(J, R) + 1j*randn(J, R));
    C = C_t + 0.05*(randn(M, R) + 1j*randn(M, R)); S = S_t + 0.05*(randn(K, R) + 1j*randn(K, R));
    for iter = 1:max_iter
        A_old = A;
        V = (S'*S) .* (C'*C) .* (B'*B); A = (Y1 * conj(khatri_rao(S, khatri_rao(C, B)))) / V.';
        V = (S'*S) .* (C'*C) .* (A'*A); B = (Y2 * conj(khatri_rao(S, khatri_rao(C, A)))) / V.';
        V = (S'*S) .* (B'*B) .* (A'*A); C = (Y3 * conj(khatri_rao(S, khatri_rao(B, A)))) / V.';
        V = (C'*C) .* (B'*B) .* (A'*A); S = (Y4 * conj(khatri_rao(C, khatri_rao(B, A)))) / V.';
        for r = 1:R, A(:,r) = A(:,r)/norm(A(:,r)); B(:,r) = B(:,r)/norm(B(:,r)); C(:,r) = C(:,r)/norm(C(:,r)); end
        if norm(A - A_old, 'fro') / norm(A, 'fro') < tol, break; end
    end
end