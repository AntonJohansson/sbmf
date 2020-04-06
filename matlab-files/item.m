t0 = 0;
tf = 500;
dt = 1;

nmeas = 1;

Lx=10*pi;
N=256;

dk = 2*pi/Lx;
dx = Lx/N;
x = -Lx/2:dx:Lx/2-dx;
k = [0:N/2-1 -N/2:-1]*2*pi/Lx;
t = t0:dt:tf;

Vx = @(x,u) 3*cos(x).^2;
Vk = @(k) 0.5 * k.^2;

Uxh = @(x,u) 	exp(-(Vx(x,u)*dt)/2);
Ux  = @(x,u) 	exp(-(Vx(x,u)*dt));
Uf  = @(k) 		exp(-(Vk(k)*dt));

phi0 = (1/Lx)*ones(size(x));
phi0 = exp(-x.^2/Lx);

phi = phi0;
phi = phi.*Uxh(x,phi);

psiF = fft(phi)*dx;
psiF = psiF.*Uf(k);
phi = ifft(psiF)/Lx;

nt = length(t);

figure(1);
for i = 1:nt
	phi = phi .* Ux(x,phi);
	psiF = fft(phi)*dx;
	psiF = psiF.*Uf(k);
	phi = ifft(psiF)/Lx;

	nn = sum(conj(phi) .* phi)*dx;
	phi = 1/sqrt(nn)*phi;

	if mod(i, nmeas) == 0
		phit = phi .* Uxh(x,phi);

		clf;
		plot(x, Vx(x,phit));
		hold on;
		plot(x, abs(phit));
		%plot(x, abs(phi_exact));
		hold off;
		title(sprintf('i: %d/%d', i,nt-1));

		drawnow;
	end
end

phi = phi .* Uxh(x,phi);
