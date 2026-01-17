function idx = find_freq_indices(f, fmin, fmax)
    idx = find(f>=fmin & f<=fmax);
    assert(~isempty(idx), 'Frequency index set is empty for [%.2f %.2f] Hz.', fmin, fmax);
end
