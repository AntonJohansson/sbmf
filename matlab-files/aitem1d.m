A=1;
P=3;
Lx=2*pi;
N=128;
c=0.7;
DT=1;
max_iteration=10000;
error_tolerance=1e-12;

dx=Lx/N;
x=-Lx/2:dx:Lx/2-dx;
kx=[0:N/2-1 -N/2:-1]*2*pi/Lx;

Uexact = pi^(-0.25) * exp(-(x).^2/2);
U=exp(-(x.^2));
U=(1/Lx)*ones(size(x));
U = U/max(max(abs(U)))*A;
%U = U/sqrt(sum(sum(abs(U).^2))*dx/P);

%potential = @(X,Y,U) 3*(cos(X).^2 + cos(Y).^2).*U + U.^3;
%potential = @(X,Y,U) 0.1*(X.^2 + Y.^2)-100;
%potential = @(X,Y) 3*(cos(X).^2 + cos(Y).^2);
%potential = @(X,Y) zeros(size(X));
Vx = 0.5*x.^2;
Vk = 0.5*kx.^2;

nmeas = 100;

figure(1);
for nn=1:max_iteration
	Uold = U;

	L00U = ifft(-(Vk).*fft(U)) + Vx.*U;

	MinvU=ifft(fft(U)./(c+Vk));
	mu=sum(sum(L00U.*MinvU))/sum(sum(U.*MinvU));

	U = U + ifft(fft(L00U-mu*U)./(c + Vk))*DT;
	U = U/max(max(abs(U)))*A;
	%U = U/sqrt(sum(sum(abs(U).^2))*dx/P);

	Uerror(nn)=sqrt(sum(sum(abs(U-Uold).^2))*dx);
	if Uerror(nn) < error_tolerance
		break
	end


	if mod(nmeas, nn) == 0
		clf;
		plot(x, Vx);
		ylim([0,1]);
		xlim([-5,5]);
		hold on;
		plot(x, abs(U));
		plot(x, abs(Uexact));
		hold off;
		drawnow;
	end
end
