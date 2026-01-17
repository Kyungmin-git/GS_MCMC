function scale = read_scale_from_meta(metaFile)
    fid = fopen(metaFile,'r');
    assert(fid>0, 'Cannot open metadata file: %s', metaFile);
    C = textscan(fid, ['%s',repmat('%s',[1,20])], 'Delimiter','|');
    fclose(fid);

    % Your original: scale = str2num(Metaseis{1,12}{3,1});
    % Add a safety check:
    raw = C{1,12}{3,1};
    scale = str2double(raw);
    assert(isfinite(scale) && scale~=0, 'Invalid scale parsed from metadata.');
end
