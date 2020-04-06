A=1.0384;
Lx=10*pi;
Ly=10*pi;
N=128;
c=0.7;
DT=1;
max_iteration=100;
error_tolerance=1e-10;

dx=Lx/N;
dy=Ly/N;
x=-Lx/2:dx:Lx/2-dx;
y=-Ly/2:dy:Ly/2-dy;
kx=[0:N/2-1 -N/2:-1]*2*pi/Lx;
ky=[0:N/2-1 -N/2:-1]*2*pi/Ly;

[X,Y]=meshgrid(x, y); [KX,KY]=meshgrid(kx, ky);

U=exp(-(X.^2 + Y.^2));
U = U/max(max(abs(U)))*A;

%potential = @(X,Y,U) 3*(cos(X).^2 + cos(Y).^2).*U + U.^3;
%potential = @(X,Y,U) 0.1*(X.^2 + Y.^2)-100;
%potential = @(X,Y) 3*(cos(X).^2 + cos(Y).^2);
%potential = @(X,Y) zeros(size(X));
potential = @(X,Y) 0.5*(X.^2 + Y.^2);
K2 = 0.5*(KX.^2 + KY.^2);

figure(1);
subplot(4,1,1);
mesh(X,Y,potential(X,Y));
for nn=1:max_iteration
	Uold = U;

	L00U = ifft2(-(K2).*fft2(U)) + potential(X,Y).*U + U.^3;

	MinvU=ifft2(fft2(U)./(c+K2));
	mu=sum(sum(L00U.*MinvU))/sum(sum(U.*MinvU));

	U = U + ifft2(fft2(L00U-mu*U)./(c + K2))*DT;
	U = U/max(max(abs(U)))*A;

	subplot(3,1,2);
	mesh(X, Y, real(U));

	Uerror(nn)=sqrt(sum(sum(abs(U-Uold).^2))*dx*dy);
	if Uerror(nn) < error_tolerance
		break
	end

	subplot(3,1,3);
	plot(1:nn, Uerror(1:nn));

	drawnow;
end
