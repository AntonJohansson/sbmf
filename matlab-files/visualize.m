N = 128;
data = dlmread('aitem_wf', '\t');
err = dlmread('aitem_error_output', '\t');

figure(1);

for i = 1:size(data,1)
	subplot(3,1,1);
	plot(1:i, err(1:i));
	title(sprintf('Error (%.3f)', err(i)));

	subplot(3,1,2);
	z = reshape(data(i,:), N, N);
	mesh(real(z));
	title("Ground state w.f. real part");

	subplot(3,1,3);
	z = reshape(data(i,:), N, N);
	mesh(imag(z));
	title("Ground state w.f. imag part");

	drawnow;
end
