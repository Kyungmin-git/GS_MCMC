function plot_complex_scatter(noiseComplexSamples, f, fTarget)
    [~,ix] = min(abs(f - fTarget));
    z = noiseComplexSamples(ix,:);
    figure; scatter(real(z), imag(z), 10, 'filled');
    axis equal; grid on;
    title(sprintf('Noise complex scatter near %.2f Hz', f(ix)));
    xlabel('Real'); ylabel('Imag'); set(gca,'FontSize',14);
end
