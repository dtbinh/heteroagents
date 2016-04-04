% Generate Rouwenhorst transition matrix for AR(1) process
function [y, Pi] = Rouwenhorst(rho, sigma, N)
	q = (1 + rho) / 2;
	p = q;
	Pi = Theta(p, q, N);
	Phi = (sigma ^ 2 / (1 - rho ^ 2)) ^ 0.5 * ((N - 1) ^ 0.5);
	y = linspace(-Phi, Phi, N)';
end

function a = Theta(p, q, n)
	if (n == 2)
		a = [[p, 1 - p]; [1 - q, q]];
    else
		b = Theta(p, q, n - 1);
		c_1 = zeros(n, n);
		c_2 = zeros(n, n);
		c_3 = zeros(n, n);
		c_4 = zeros(n, n);
		c_1(1:n - 1, 1:n - 1) = b;
		c_2(1:n - 1, 2:n) = b;
		c_3(2:n, 1:n - 1) = b;
		c_4(2:n, 2:n) = b;
		a = p * c_1 + (1 - p) * c_2 + (1 - q) * c_3 + q * c_4;
		a(2:n - 1, :) = a(2:n - 1, :) / 2;
    end
end
