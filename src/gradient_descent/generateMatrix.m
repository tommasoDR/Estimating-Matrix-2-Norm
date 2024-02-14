% Function that generate a m*n matrix with desired condition number.
%
% Input:
% - m: first dimension of the matrix
% - n: second dimension of the matrix
% - upperlimit: value of the biggest element in the matrix
% - desiredCondition: desired condition number of the matrix
%
% Output:
% - A: the m*n matrix generated

function [A] = generateMatrix(m, n, upperlimit, desiredCondition)
    % Generate a random orthogonal m*m matrix U
    [U, ~] = qr(rand(m, m));

    % Generate a random orthogonal n*n matrix V
    [V, ~] = qr(rand(n, n));

    % Generate the singular values
    singularValues = linspace(desiredCondition, 1, min(m,n));

    % Create the diagonal matrix S with the singular values
    S = zeros(m,n);
    S(1:min(m,n),1:min(m,n)) = diag(singularValues);

    % Calculate the matrix A using SVD
    A = U * S * V';

    % Normalize the matrix A_svd to achieve the desired condition number
    A = A * (upperlimit/max(max(A)));
end

