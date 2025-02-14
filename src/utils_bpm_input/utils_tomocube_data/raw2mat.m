function [data] = raw2mat(infile, arrsize)
    mhdpath = strrep(infile, 'raw', 'mhd');
    if exist(mhdpath)
        arrsize = MhdDimension(mhdpath);
    end
    
    elements = prod(arrsize);

    s = dir(infile);
    bytes = s.bytes;
    
    data = [];
    type = '';
    if bytes == elements
        type = 'uint8';
    else if bytes == (2*elements)
         type = 'uint16';
    else
        type = 'single';
    end
    
    fid = fopen(infile);
    data = fread(fid, elements, type);
    fclose(fid);

    data = reshape(data, arrsize);
end
